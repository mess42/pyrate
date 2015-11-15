#!/usr/bin/env/python
"""
Pyrate - Optical raytracing based on Python

Copyright (C) 2014 Moritz Esslinger moritz.esslinger@web.de
               and Johannes Hartung j.hartung@gmx.net
               and    Uwe Lippmann  uwe.lippmann@web.de

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

class OptimizableVariable(object):
    """
    Class that contains a float. Used to get a pointer on a variable.
    """
    def __init__(self, variable_status=False, variable_type="Variable", **kwargs):
        """
        Name is gone since it's only needed for reference. Therefore the former listOfOptimizableVariables
        will be a dictionary. type may contain a string e.g.
        "Variable", "Solve", "Pickup", "...",
        status=True means updating during optimization run
        status=False means no updating during optimization run
        kwargs depend on type
        "Variable" is value=value
        "Solve" gets function=f, args=tupleof optimizablevariables
        """
        # TODO: pickup gets function=f, args=tupleofexternal variables
        self.__var_type = variable_type.lower()
        self.status = variable_status
        self.parameters = kwargs

    @property
    def var_type(self):
        return self.__var_type.lower()

    @var_type.setter
    def var_type(self, variable_type="Variable"):
        self.__var_type = variable_type.lower()

    def setvalue(self, value):
        if self.var_type == "variable":
            self.parameters["value"] = value

    def evaluate(self):
        # notice: evaluation code is not limited to floats
        if self.var_type == "variable":
            # if type = variable then give only access to value
            return self.parameters["value"]
        elif self.var_type == "solve":
            # if type = variable then pack all arguments into one tuple
            # and put it into the userdefined function
            # evaluate the result
            arguments_for_function_eval = (argfunc.evaluate() for argfunc in self.parameters["args"])
            return self.parameters["function"](*arguments_for_function_eval)

class ClassWithOptimizableVariables(object):
    """
    Implementation with dictionary.
    """
    def __init__(self):
        self.dict_variables = {}

    def addVariable(self, name, var):
        self.dict_variables[name] = var

    def getAllVariables(self):
        """
        Conversion of dict into list of Variables. These are only references to the objects in dict.
        Therefore all changes to them directly affects the dictionary in the class.
        """
        # TODO: if there are sub classes which are inherited from ClassWith... append their getAllVariables()
        # TODO: since these are also containing references to the objects in their respective dicts there should
        # TODO: not be any problem that there are copies created

        lst_of_vars = self.dict_variables.values()
        lst_of_attributes_which_are_class_with_opt_vars = filter(lambda x: isinstance(x, ClassWithOptimizableVariables), self.__dict__.values())

        for a in lst_of_attributes_which_are_class_with_opt_vars:
            lst_of_vars.extend(a.getAllVariables())

        return lst_of_vars

    def getAllValues(self):
        """
        For fast evaluation of value vector
        """
        return np.array([a.evaluate() for a in self.getAllVariables()])

    def getAllStates(self):
        """
        For fast evaluation of optimizable values vector
        """
        return np.array([a.status for a in self.getAllVariables()])

    def getActiveVariables(self):
        """
        Returns a list of active variables the names are lost
        but it does not matter since the variable references are still in the
        original dictionary in the class
        """
        return filter(lambda x: x.status, self.getAllVariables())

    def getActiveValues(self):
        return np.array([a.evaluate() for a in self.getActiveVariables()])

    def setStatus(self, name, var_status=True):
        self.dict_variables[name].status = var_status

class ExampleSubClass(ClassWithOptimizableVariables):
    def __init__(self):
        self.a = 10.

class ExampleSuperClass(ClassWithOptimizableVariables):
    def __init__(self):
        super(ExampleSuperClass, self).__init__()
        self.b = 20.
        self.c = ClassWithOptimizableVariables()
        self.c.addVariable("blubberbla", OptimizableVariable(False, "Variable", value=5.0))
        self.addVariable("blubberdieblub", OptimizableVariable(False, "Variable", value=10.0))




if __name__ == "__main__":
    def f(p, q):
        return p + q

    p = OptimizableVariable(False, "Variable", value="glass1")
    q = OptimizableVariable(True, "Variable", value="glass2")
    r = OptimizableVariable(False, "Solve", function=f, args=(p, q))
    print p.__dict__
    print p.evaluate()
    print q.__dict__
    print q.evaluate()
    print r.__dict__
    print r.evaluate()
    p.setvalue("glass5") # TODO: assignment operator overloading
    q.setvalue("glass6") # TODO: should behave different for Variable or Solve
    print r.evaluate()

    cl = ClassWithOptimizableVariables()
    cl.addVariable("var1", p)
    cl.addVariable("var2", q)
    cl.addVariable("var3", r)
    print cl.__dict__
    print cl.getAllVariables()
    print cl.getAllValues()
    print cl.getAllStates()
    lst = cl.getActiveVariables()
    print lst[0].evaluate()
    lst[0].setvalue("blublub")
    print cl.getAllValues()[cl.getAllStates()]
    print cl.getAllValues()
    print cl.getActiveValues()

    cl2 = ExampleSuperClass()
    print cl2.__dict__
    print cl2.getAllVariables()
    print cl2.getAllValues()



