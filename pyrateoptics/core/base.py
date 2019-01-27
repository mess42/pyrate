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

import numpy as np
import math
import ruamel.yaml as yaml
import re

from copy import copy

from .log import BaseLogger


class OptimizableVariable(BaseLogger):
    """
    Class that contains an optimizable variable.
    Used to get a pointer on a variable.
    The value is not constrained to float.
    Also other dependent variables are possible to define.
    """
    def __init__(self, variable_type="fixed",
                 name="", kind="optimizablevariable", **kwargs):

        super(OptimizableVariable, self).__init__(
                name=name,
                kind=kind,
                **kwargs)

        """
        kwargs depend on type
        "variable" is value=value
        "pickup" gets function=f, args=tuple of optimizablevariables
        "external" gets function=f, args=tuple of external standard variables
        (or values)

        :param variable_type (string): which kind of variable do we have?
        valid choices are: 'fixed', 'variable', 'pickup', 'external'

        :param kwargs (keyword arguments): used to fill up parameters
        dictionary for variable.

        kwargs:

        type:        kwargs:
        'fixed'         value=1.0 or value='hello'
        'variable'    value=1.0, or value='hello'
        'pickup'        function=f, args=(other_optvar1, other_optvar2, ...)
        'external'        function=f, args=(other_externalvar/value, ...)

        Notice: f must take exactly as many arguments as args=... long is.
        The only constraint on f is that it has to convert the variables into
        some final float variable. If this is not the case some
        low-level optimizer might break down.


        """

        self.evaldict = {
                    "fixed": self.eval_fixed,
                    "variable": self.eval_variable,
                    "pickup": self.eval_pickup,
                    "external": self.eval_external
                    }

        self.initdict = {
                    "fixed": self.init_fixed,
                    "variable": self.init_variable,
                    "pickup": self.init_pickup,
                    "external": self.init_external
                    }

        self.var_type = variable_type
        self.evalfunc = self.evaldict[self.var_type]
        self.initdict[self.var_type](**kwargs)
        self.set_interval(None, None)

    def init_fixed(self, **kwargs):
        self.parameters = {}
        self.parameters["value"] = kwargs.get("value", None)

    def init_variable(self, **kwargs):
        self.parameters = {}
        self.parameters["value"] = kwargs.get("value", None)

    def init_pickup(self, **kwargs):
        self.parameters = {}
        self.parameters["functionobject"] = kwargs.get("functionobject",
                                                       (None, ()))
        # self.parameters["function"] = kwargs.get("function", None)
        self.parameters["args"] = kwargs.get("args", ())
        # TODO: function also as string, string tuple or as code object

    def init_external(self, **kwargs):
        self.parameters = {}
        self.parameters["functionobject"] = kwargs.get("functionobject",
                                                       (None, ()))
        self.parameters["args"] = kwargs.get("args", ())
        # TODO: function also as string, string tuple or as code object

    def getVariableType(self):
        return self.__var_type.lower()

    def setVariableType(self, vtype):
        self.__var_type = vtype.lower()

    var_type = property(fget=getVariableType, fset=setVariableType)

    def set_interval(self, left=None, right=None):
        """
        Interval transform from finite to infinite and back.
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
        # TODO: fine-tune for one-sided intervals
        if left is None and right is None:
            self.transform = lambda x: x
            self.inv_transform = lambda x: x
        else:
            self.inv_transform = lambda x: left +\
                (right - left)/(1. + math.exp(-x/math.fabs(right - left)))
            self.transform = lambda x: math.log((-x + left)/(x - right)) *\
                math.fabs(left - right)

    """
    all -> fixed: value is conserved
    all -> variable: value is conserved
    all -> pickup: value is not conserved
    all -> external: value is not conserved
    """

    def changetype_from_to(self, from_type, to_type, **kwargs):
        """
        Function only called indirectly from changetype() to make program
        logic more clear.
        """

        # first: backup value
        # second: backup parameters
        # third: erase parameters (will be done in init functions)
        # reset parameters: also done in init functions

        self.debug("changing type from \'%s\' to \'%s\'" % (from_type,
                                                            to_type))
        value_backup = self.evaluate()
        parameters_backup = self.parameters
        self.debug("old value (%s) and old parameters (%s)" %
                   (str(value_backup), str(parameters_backup)))
        self.debug("kwargs=%s" % repr(kwargs))
        self.initdict[to_type](**kwargs)

        if to_type == "fixed" or to_type == "variable":
            if "value" not in kwargs:
                self.parameters["value"] = value_backup

        if to_type == "pickup" or to_type == "external":
            (functionobject, functionname) = self.parameters["functionobject"]
            functionobject.generateFunctionsFromSource([functionname])

        if (from_type == "pickup" and to_type == "external") or\
                (from_type == "external" and to_type == "pickup"):
            #  if not kwargs.has_key("function"):
            #    self.parameters["function"] = parameters_backup["function"]
            if "functionobject" not in kwargs:
                self.parameters["functionobject"] =\
                    parameters_backup["functionobject"]
            self.info("PICKUP SET")

        self.evalfunc = self.evaldict[to_type]
        self.var_type = to_type
        self.debug("new value %s and new parameters %s" %
                   (str(self.evaluate()), str(self.parameters)))

    def changetype(self, vtype, **kwargs):
        self.changetype_from_to(self.var_type, vtype.lower(), **kwargs)

    def setvalue(self, value):
        # TODO: overload assign operator
        if self.var_type == "variable" or self.var_type == "fixed":
            self.parameters["value"] = value

    def eval_fixed(self):
        # if type = variable then give only access to value
        return self.parameters.get("value", None)

    def eval_variable(self):
        # if type = variable then give only access to value
        return self.parameters.get("value", None)

    def eval_pickup(self):
        # if type = pickup then package up all arguments into one tuple
        # and put it into the userdefined function
        # evaluate the result
        arguments_for_function_eval = (argfunc.evaluate()
                                       for argfunc in self.parameters["args"])
        (functionobject, functionname) = self.parameters["functionobject"]
        return functionobject.functions[functionname](
                *arguments_for_function_eval)

    def eval_external(self):
        # same as for solve except that there are no further
        # OptimizableVariables to be considered
        (functionobject, functionname) = self.parameters["functionobject"]
        return functionobject.functions[functionname](*self.parameters["args"])

    def evaluate(self):
        # notice: evaluation code is not limited to floats
        # do not use if-else construct because for many
        # cases dict lookup is faster

        return self.evalfunc()

    def __call__(self):
        return self.evaluate()

    def evaluate_transformed(self):
        """
        Transform variable value from finite interval
        to infinite IR before evaluation
        """
        return self.transform(self.evaluate())

    def setvalue_transformed(self, value_transformed):
        """
        Transform value back from infinite IR before setting value
        """
        value = self.inv_transform(value_transformed)
        self.setvalue(value)

    def getDictionary(self):
        res = super(OptimizableVariable, self).getDictionary()
        res["variable_type"] = self.var_type
        for (key, val) in self.parameters.items():
            if self.var_type == 'pickup':
                if key == 'args':
                    val = tuple([v.unique_id for v in val])
                if key == 'functionobject':
                    val = (val[0].source, val[1])
            elif self.var_type == "fixed" or self.var_type == "variable":
                if key == "value":
                    val = float(val)
                # TODO: what about other types?
            res[key] = val
        return res

    @staticmethod
    def initFromDictionary(override_unique_id=True, **kwargs):
        """
        Creates OptimizableVariable from dictionary entry from serialization
        dictionary. Overriding the unique_id is standard behaviour since
        unique_id is only needed for link purposes and is not intended to be
        stable.
        """
        res = OptimizableVariable(**kwargs)
        if not override_unique_id:
            res.unique_id = kwargs["unique_id"]
        return res

# TODO: This class contains far too much stuff! Refactor!
# Four core functionalities which may be splitted:
#   I   Observer (lightweight) - needed for interface communication
#   II  Collecting sub classes and variables (heavy) - needed for III and IV
#   III Provide an API for preparing a numpy array for the optimizer (heavy)
#   IV  Some serialization techniques


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
    def __init__(self, name="", kind="classwithoptimizablevariables",
                 **kwargs):
        """
        Initialize with empty dict.
        """
        super(ClassWithOptimizableVariables, self).__init__(
                name=name,
                kind=kind,
                **kwargs)

        self.annotations = {}
        self.list_observers = []
        # for the optimizable variable class it is useful to have some observer
        # links they get informed if variables change their values

    def appendObservers(self, obslist):
        self.list_observers += obslist

    def informObservers(self):
        for obs in self.list_observers:
            obs.informAboutUpdate()

    def getDictionary(self):
        res = super(ClassWithOptimizableVariables, self).getDictionary()
        res["annotations"] = self.annotations
        res["variables"] = self.getTypesForDict(typ=OptimizableVariable)
        res["classes"] = self.getTypesForDict(typ=ClassWithOptimizableVariables)
        return res

    def getDictionaryAllVariablesById(self):
        return dict([(v.unique_id, v.getDictionary()) for v in self.getAllVariables()["vars"].values()])

    def getDictionaryAllClassesById(self, removeself=True):
        res = dict([(v.unique_id, v.getDictionary()) for v in self.getAllVariables()["classes"].values()])
        if removeself:
            res.pop(self.unique_id)
        return res

    def getDictionaryClassesById(self):
        return dict([(v.unique_id, v.getDictionary())
                     for v in self.getTypesForDict(ClassWithOptimizableVariables, lambda x: x).values()])

    def getTypesForDict(self, typ, func=lambda x: x.unique_id):

        def remove_non_optvars(var):
            if isinstance(var, typ):
                    return func(var)
            elif isinstance(var, list):
                l1 = []
                for item in var:
                    newitem = remove_non_optvars(item)
                    if newitem is not None:
                        l1.append(newitem)
                if l1 != []:
                    return l1
                else:
                    return None
            elif isinstance(var, dict):
                d1 = {}
                for (key, value) in var.items():
                    newitem = remove_non_optvars(value)
                    if newitem is not None:
                        d1[key] = newitem
                if d1 != {}:
                    return d1
                else:
                    return None

            return None

        mydict = {}

        for (key, value) in self.__dict__.items():
            myitem = remove_non_optvars(value)
            if myitem is not None:
                mydict[key] = myitem

        return mydict

    def getAllVariables(self, recursive=True):

        """
        This functions traverses through a ClassWithOptimizableVariables
        and collects several information of the underlying recursive structure.
        """

        # TODO: Rename function to a more appropriate name. (form vs. function)
        # TODO: This function is too intelligent. Divide into more appropriate sub problems

        def addOptimizableVariablesToList(var,
                                          dictOfOptVars={"vars": {},
                                                         "longkeystrings": {},
                                                         "deref": {},
                                                         "classes": {}},
                                          idlist=[],
                                          keystring="",
                                          reducedkeystring="", recursive=True):
            """
            Accumulates optimizable variables in var and its linked objects.
            Ignores ring-links and double links.

            Variables are accumulated under the "vars" key.
            Longkeystrings are accumulated under the "longkeystrings" key.
            Strings for object dereference by eval are accumulated under the
            "deref" key.

            Below those main keys there is a dict in which the key is a
            combination of names of the variables for an
            easy-to-remember-reference to the variables.

            @param var: object to evaluate (object)
            @param dictOfOptVars: optimizable variables found (dict of dict of objects)
            @param idlist: ids of objects already evaluated (list of int)
            """

            if id(var) not in idlist:
                idlist.append(id(var))

                if isinstance(var, ClassWithOptimizableVariables):
                    if len(idlist) == 1 or (len(idlist) > 1 and recursive):
                        dictOfOptVars["classes"][var.unique_id] = var
                        for (k, v) in var.__dict__.items():
                            newkeystring = keystring + "(" + var.name + ")." +\
                                str(k)
                            newredkeystring = reducedkeystring + var.name + "."
                            dictOfOptVars, idlist =\
                                addOptimizableVariablesToList(v,
                                                              dictOfOptVars,
                                                              idlist,
                                                              newkeystring,
                                                              newredkeystring,
                                                              recursive=recursive)
                elif isinstance(var, dict):
                    for (k, v) in var.items():
                        newkeystring = keystring + "[\"" + str(k) + "\"]"
                        dictOfOptVars, idlist =\
                            addOptimizableVariablesToList(v,
                                                          dictOfOptVars,
                                                          idlist,
                                                          newkeystring,
                                                          reducedkeystring,
                                                          recursive=recursive)

                elif isinstance(var, list) or isinstance(var, tuple):
                    for (ind, v) in enumerate(var):
                        newkeystring = keystring + "[" + str(ind) + "]"
                        dictOfOptVars, idlist =\
                            addOptimizableVariablesToList(v,
                                                          dictOfOptVars,
                                                          idlist,
                                                          newkeystring,
                                                          reducedkeystring,
                                                          recursive=recursive)

                elif isinstance(var, OptimizableVariable):
                    newkeystring = keystring + "(" + var.name + ")"
                    newreducedkeystring = reducedkeystring + var.name
                    dictOfOptVars["vars"][newreducedkeystring] = var
                    dictOfOptVars["longkeystrings"][newreducedkeystring] =\
                        newkeystring
                    derefstring = re.sub(r"\([a-zA-Z0-9_ ]+\)", "",
                                         newkeystring)
                    dictOfOptVars["deref"][newreducedkeystring] = derefstring

            return dictOfOptVars, idlist

        """
        Traverse ClassWithOptimizableVariables and their subclasses
        and substructures. It does not work for OptimizableVariables which
        are arguments to a pickup and were defined outside of the actual
        ClassWithOptimizableVariables.
        """

        (dict_opt_vars, idlist) =\
            addOptimizableVariablesToList(self,
                                          keystring="",
                                          recursive=recursive)

        """
        Fix the former neglection of pickup variables from outside self.
        """
        variables_without_pickup_arguments = [(k, v) for (k, v) in dict_opt_vars["vars"].items()]

        for (key, var) in variables_without_pickup_arguments:
            if var.var_type == "pickup":
                for (ind, v) in enumerate(var.parameters.get("args", [])):
                    newkeystring = key + "[" + str(ind) + "]"
                    newreducedkeystring = key + var.name + ".pickup_args[" + str(ind) + "]="
                    (dict_opt_vars, idlist) =\
                        addOptimizableVariablesToList(v,
                                                      dict_opt_vars,
                                                      idlist,
                                                      newkeystring,
                                                      newreducedkeystring,
                                                      recursive=recursive)




        return dict_opt_vars

        # TODO: due to the dict the order of the variables is sometimes not
        # maintained between to consecutive runs of the programs; this is
        # maybe not good

    def getAllValues(self):
        """
        For fast evaluation of value vector
        """
        return np.array([a.evaluate() for a in self.getAllVariables()["vars"].
                         values()])

    def getActiveVariables(self):
        """
        Returns a list of active variables the names are lost
        but it does not matter since the variable references are still in the
        original dictionary in the class
        """
        return [x for x in self.getAllVariables()["vars"].values()
                if x.var_type == "variable"]

    def getActiveValues(self):
        """
        Function to get all values into one large np.array.
        """
        return np.array([a() for a in self.getActiveVariables()])

    def setActiveValues(self, x):
        """
        Function to set all values of active variables to the values in the
        large np.array x.
        """
        for i, var in enumerate(self.getActiveVariables()):
            var.setvalue(x[i])

    def getActiveTransformedValues(self):
        return np.array([a.evaluate_transformed() for a in
                         self.getActiveVariables()])

    def setActiveTransformedValues(self, x):
        """
        Function to set all values of active variables to the values in the
        large np.array x.
        """
        for i, var in enumerate(self.getActiveVariables()):
            var.setvalue_transformed(x[i])

    def resetVariable(self, key, var):
        """
        Resets variable by evaluating deref string from getAllVariables
        dictionary. Maybe this could be solved by using less black magic.
        Although this method is necessary for maintaining multi configs.
        """

        dict_of_vars = self.getAllVariables()
        deref = dict_of_vars["deref"][key]
        exec(compile("self" + deref + " = var", "<string>", "exec"))

    def getVariable(self, key):
        """
        Gets variable from short key.
        """
        dict_of_vars = self.getAllVariables()
        variable = dict_of_vars["vars"][key]
        return variable

    def getCompleteListForReconstruction(self):
        return [self.getDictionary(),
                self.getDictionaryAllClassesById(),
                self.getDictionaryAllVariablesById()]
