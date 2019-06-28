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


from .functionobject import FunctionObject
from .log import BaseLogger


class State(object):
    """
    State basic object which uses the State pattern of the gang of four.
    It allows to control transitions between states by implementing them
    for every State derived class and implementing an interface by using
    a single object. This is to keep large transition if-then-else
    constructions minimal.
    """
    def __init__(self):
        self.var_type = "abstract"
        self.parameters = {}

    def evaluate(self):
        raise NotImplementedError()

    def __call__(self):
        return self.evaluate()

    def set_value(self, value):
        raise NotImplementedError()

    def toFixed(self, context):
        raise NotImplementedError()

    def toVariable(self, context):
        raise NotImplementedError()

    def toPickup(self, context, functionobject_functionname_tuple, args):
        raise NotImplementedError()

    def toDictionary(self):
        raise NotImplementedError()

    @staticmethod
    def fromDictionary(dictionary, *args):
        raise NotImplementedError()


class SimpleValueState(State):
    """
    Base class for FixedState and VariableState since they only
    differ by type.
    """
    def __init__(self, var_type, value):
        super(SimpleValueState, self).__init__()
        self.var_type = var_type
        self.parameters["value"] = value

    def set_value(self, value):
        self.parameters["value"] = value

    def evaluate(self):
        return self.parameters["value"]

    def toFixed(self, context):
        context._state = FixedState(self.parameters["value"])

    def toVariable(self, context):
        context._state = VariableState(self.parameters["value"])

    def toPickup(self, context, functionobject_functionname_tuple, args):
        context._state = PickupState(functionobject_functionname_tuple, args)

    def toDictionary(self):
        resdict = {}
        resdict["value"] = self.parameters["value"]
        return resdict


class FixedState(SimpleValueState):
    def __init__(self, value):
        super(FixedState, self).__init__("fixed", value)

    @staticmethod
    def fromDictionary(dictionary, *args):
        return FixedState(dictionary.pop("value"))


class VariableState(SimpleValueState):
    def __init__(self, value):
        super(VariableState, self).__init__("variable", value)

    @staticmethod
    def fromDictionary(dictionary, *args):
        return VariableState(dictionary.pop("value"))


class PickupState(State):
    def __init__(self, functionobject_functionname_tuple, args, isvalid=True):
        """
        @param functionobject_functionname_tuple
            (functionobject,
             functionname)

        @param args (list, tuple) of other OptimizableVariables

        @param isvalid (boolean) whether initialization was completed

        Notice: function "functionname" must take exactly as many arguments
        as args long is.
        The only constraint on f is that it has to convert the variables into
        some final variable. If this is not the case some
        low-level optimizer might break down.
        """

        super(PickupState, self).__init__()
        self.var_type = "pickup"
        (functionobject, functionname) = self.parameters["functionobject"] =\
            functionobject_functionname_tuple
        self.parameters["args"] = args
        self.isvalid = isvalid
        functionobject.generateFunctionsFromSource([functionname])

    def set_value(self, value):
        pass

    def toFixed(self, context):
        context._state = FixedState(self.evaluate())

    def toVariable(self, context):
        context._state = VariableState(self.evaluate())

    def toPickup(self, context, functionobject_functionname_tuple, args):
        context._state = PickupState(functionobject_functionname_tuple, args)

    def evaluate(self):
        if not self.isvalid:
            return None
        arguments_for_function_eval = (argfunc.evaluate()
                                       for argfunc in self.parameters["args"])
        (functionobject, functionname) = self.parameters["functionobject"]
        return functionobject.functions[functionname](
                *arguments_for_function_eval)

    def toDictionary(self):
        resdict = {}
        (functionobject, functionname) = self.parameters["functionobject"]
        resdict["functionobject"] = functionobject.toDictionary()
        resdict["functionname"] = functionname
        resdict["args"] = list([a.unique_id for a in self.parameters["args"]])
        return resdict

    @staticmethod
    def fromDictionary(dictionary, *args):
        (source_checked, variables_checked) = args
        fo = FunctionObject.fromDictionary(dictionary["functionobject"],
                                           source_checked, variables_checked)
        functionname = dictionary["functionname"]
        # isvalid=False is set because args are not initialized yet
        # this is done later in the pool initialization
        return PickupState((fo, functionname), dictionary["args"],
                           isvalid=False)


class OptimizableVariable(BaseLogger):
    """
    Optimizable variable which provides an interface to the
    State objects given above. The value is not constrained to float.
    It provides also a possibility to do one inline transform.
    Also other dependent variables are possible to define.
    """

    id_trafo = (FunctionObject("identity = lambda x: x", ["identity"],
                               name="id_trafo"),
                "identity", "identity")

    def __init__(self, state, name="", unique_id=None):
        """
        @param state (State): which kind of variable do we have?
        valid choices are: FixedState, VariableState, PickupState

        kwargs:

        type:        kwargs:
        'fixed'         value=1.0 or value='hello'
        'variable'    value=1.0, or value='hello'
        'pickup'        function=f, args=(other_optvar1, other_optvar2, ...)
        'external'        function=f, args=(other_externalvar/value, ...)
        """


        super(OptimizableVariable, self).__init__(name=name,
                                                  kind="optimizablevariable",
                                                  unique_id=unique_id)
        self._state = state
        self.set_transform(self.id_trafo)

    def toFixed(self):
        """
        State transition to fixed.
        """
        self._state.toFixed(self)

    def toVariable(self):
        """
        State transition to variable.
        """
        self._state.toVariable(self)

    def toPickup(self, functionobject_functionname_tuple, args):
        """
        State transition to pickup.

        @param functionobject_functionname_tuple is tuple of
                functionobject and functionname.
        @args list or tuple of other optimizable variables
                which fit as arguments list into functionname.
        """
        self._state.toPickup(self, functionobject_functionname_tuple, args)

    def evaluate(self):
        """
        Evaluation of value of optimizable variable.
        """
        return self._state.evaluate()

    def set_value(self, value):
        """
        Set value of optimizable variable: This makes only sense for
        fixed and variable states.
        """
        self._state.set_value(value)

    def __call__(self):
        """
        Short form of evaluate.
        """
        return self.evaluate()

    def var_type(self):
        """
        Returns variable type as string
        """
        return self._state.var_type

    def set_transform(self, functionobject_functionname_triple):
        """
        Set transform together with inverse transform from one
        function object:

        @param functionobject_functionname_triple
                (functionobject, transformname, invtransformname)
        """
        self._transform_functionobject_functionname_triple =\
            functionobject_functionname_triple
        (fo, ftrafo, finvtrafo) =\
            self._transform_functionobject_functionname_triple
        fo.generateFunctionsFromSource([ftrafo, finvtrafo])

    def evaluate_transformed(self):
        """
        Evaluates to transformed value.
        """
        (fo, ft, fit) = self._transform_functionobject_functionname_triple
        return fo.functions[ft](self.evaluate())

    def set_value_transformed(self, value_transformed):
        """
        Sets transformed value by performing the inverse transform
        on value_transformed and setting the value from this.
        """
        (fo, ft, fit) = self._transform_functionobject_functionname_triple
        value = fo.functions[fit](value_transformed)
        self.set_value(value)

    def toDictionary(self):
        """
        Providing a dictionary from the OptimizableVariable
        which can be used for reconstruction.
        """
        resdict = self.getBasicInfo()
        resdict["variable_type"] = self.var_type()
        (fo, ftrafo, finvtrafo) =\
            self._transform_functionobject_functionname_triple
        resdict["transform"] = fo.toDictionary()
        resdict["transform_name"] = ftrafo
        resdict["invtransform_name"] = finvtrafo
        resdict["state"] = self._state.toDictionary()
        return resdict

    @staticmethod
    def fromDictionary(dictionary, source_checked, variables_checked):
        """
        Providing a function to reconstruct an OptimizableVariable
        from some dictionary. The variables_pool is necessary if
        pickup variables are arising.
        """
        reconstruct_dictionary = {"fixed": FixedState,
                                  "variable": VariableState,
                                  "pickup": PickupState}
        state_dictionary = dictionary["state"]
        var_type = dictionary["variable_type"]
        state = reconstruct_dictionary[var_type].fromDictionary(
                state_dictionary,
                source_checked, variables_checked)
        ov = OptimizableVariable(state, name=dictionary["name"],
                                 unique_id=dictionary["unique_id"])
        fo = FunctionObject.fromDictionary(dictionary["transform"],
                                           source_checked, variables_checked)
        ft = dictionary["transform_name"]
        fi = dictionary["invtransform_name"]
        ov.set_transform((fo, ft, fi))
        return ov


class FloatOptimizableVariable(OptimizableVariable):
    """
    Provides an interface to a definitely float variable which
    can also be mapped from a finite interval to an infinite interval
    which is useful for optimization later.

    The reason for this is that the optimizer backend has
    only to deal with infinite range variables and the user
    can transparently change the interval. The idea and
    verification is (up to a slight change in the
    transformation) taken from:
    M. Roeber, "Multikriterielle Optimierungsverfahren fuer
    rechenzeitintensive technische Aufgabenstellungen",
    Diploma Thesis, Technische Universitaet Chemnitz,
    2010-03-31
    """

    interval_trafo_source = """

import math

def left_bounded(x):
    return math.fabs(left)*math.log((x - left)/math.fabs(left))

def left_bounded_inv(x):
    return left + math.fabs(left)*math.exp(x/math.fabs(left))

def right_bounded(x):
    return -math.fabs(right)*math.log((right - x)/math.fabs(right))

def right_bounded_inv(x):
    return right - math.fabs(right)*math.exp(-x/math.fabs(right))

def both_bounded(x):
    return math.log((-x + left)/(x - right)) * math.fabs(left - right)

def both_bounded_inv(x):
    return left +\
            (right - left)/(1. + math.exp(-x/math.fabs(right - left)))

    """

    def setInterval(self, left=None, right=None):
        self.interval_trafo_fo = FunctionObject(
                initial_sourcecode=self.interval_trafo_source,
                initial_globals_dictionary={
                                            "left": left,
                                            "right": right},
                initial_function_names=["left_bounded", "left_bounded_inv",
                                        "right_bounded", "right_bounded_inv",
                                        "both_bounded", "both_bounded_inv"],
                name="interval_trafo")

        if left is not None and right is None:
            self.set_transform((self.interval_trafo_fo,
                                "left_bounded",
                                "left_bounded_inv"))
        elif left is None and right is not None:
            self.set_transform((self.interval_trafo_fo,
                                "right_bounded",
                                "right_bounded_inv"))
        elif left is not None and right is not None:
            self.set_transform((self.interval_trafo_fo,
                                "both_bounded",
                                "both_bounded_inv"))
