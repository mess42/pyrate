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
from scipy.optimize import minimize


class OptimizableVariable(object):
    """
    Class that contains an optimizable variable. Used to get a pointer on a variable.
    The value is not constrained to float. Also other dependent variables are possible to define.
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
        "Solve" gets function=f, args=tuple of optimizablevariables
        "Pickup" gets function=f, args=tuple of external standard variables (or values)

        :param variable_status (bool): should variable optimized during optimization run?
        :param variable_type (string): which kind of variable do we have? valid choices are: 'variable', 'solve', 'pickup'

        :param kwargs (keyword arguments): used to fill up parameters dictionary for variable.

        kwargs:

        type:        kwargs:
        'variable'    value=1.0, or value='hello'
        'solve'        function=f, args=(other_optvar1, other_optvar2, ...)
        'pickup'        function=f, args=(other_externalvar/value, ...)

        Notice: f must take exactly as many arguments as args=... long is. The only constraint on f is
        that it has to convert the variables into some final float variable. If this is not the case some
        low-level optimizer might break down.


        """
        self.__var_type = variable_type.lower()
        self.status = variable_status
        self.parameters = kwargs

        # TODO: status variable for "solve" and "pickup"
        # TODO: either status variable is important or not
        # TODO: most simple solution -> ignore status variable for "solve" and "pickup"
        # TODO: more complex solution -> use status variable to update initial_value and return its value
        # in evaluate()

        self.initial_value = self.evaluate()

    def __str__(self, *args, **kwargs):
        return "OptVar('" + self.var_type + "') = " + str(self.parameters)

    @property
    def var_type(self):
        return self.__var_type.lower()

    @var_type.setter
    def var_type(self, variable_type="Variable"):
        self.__var_type = variable_type.lower()

    def setvalue(self, value):
        if self.var_type == "variable":
            self.parameters["value"] = value

    def eval_variable(self):
        # if type = variable then give only access to value
        return self.parameters["value"]

    def eval_solve(self):
        # if type = variable then pack all arguments into one tuple
        # and put it into the userdefined function
        # evaluate the result
        arguments_for_function_eval = (argfunc.evaluate() for argfunc in self.parameters["args"])
        return self.parameters["function"](*arguments_for_function_eval)

    def eval_pickup(self):
        # same as for solve except that there are no further OptimizableVariables to be considered
        return self.parameters["function"](*self.parameters["args"])



    def evaluate(self):
        # notice: evaluation code is not limited to floats
        # do not use if-else construct because for many cases dict lookup is faster
        evaldict = {
                    "variable" : self.eval_variable,
                    "solve" : self.eval_solve,
                    "pickup" : self.eval_pickup
                    }

        return evaldict[self.var_type]()

class ClassWithOptimizableVariables(object):
    """
    Implementation of some class with optimizable variables with the help of a dictionary.
    This class is also able to collect the variables and their values from its subclasses per recursion.
    """
    def __init__(self):
        """
        Initialize with empty dict.
        """
        self.dict_variables = {}

    def addVariable(self, name, var):
        """
        Add some variable into dict.
        """
        self.dict_variables[name] = var

    def getAllVariables(self):
        """
        Conversion of dict into list of Variables. These are only references to the objects in dict.
        Therefore all changes to them directly affects the dictionary in the class.
        """
        # if there are sub classes which are inherited from ClassWith... append their getAllVariables()
        # since these are also containing references to the objects in their respective dicts there should
        # not be any problem that there are copies created

        lst_of_vars = self.dict_variables.values()
        lst_of_attributes_which_are_class_with_opt_vars = filter(lambda x: isinstance(x, ClassWithOptimizableVariables), self.__dict__.values())
        for list_vars in filter(lambda x: isinstance(x, list) or isinstance(x, tuple), self.__dict__.values()):
            for a in list_vars:
                lst_of_vars.extend(a.getAllVariables())

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

    def getActiveVariablesIsWriteable(self):
        """
        Returns a np.array of bools which is True for a variable and False for a solve or pickup,
        since they cannot be modified. The length of this vector is the same as for getActiveVariables()
        """
        return np.array(filter(lambda x: x.var_type == "variable", self.getActiveVariables()))


    def getActiveValues(self):
        """
        Function to get all values into one large np.array.
        """
        return np.array([a.evaluate() for a in self.getActiveVariables()])

    def setActiveValues(self, x):
        """
        Function to set all values of active variables to the values in the large np.array x.
        """
        for i, var in enumerate(self.getActiveVariables()):
            var.setvalue(x[i])

    def getActiveVariableValues(self):
        """
        Function to get all values from active variables into one large np.array.
        """
        return np.array([a.evaluate() for a in self.getActiveVariablesIsWriteable()])

    def setActiveVariableValues(self, x):
        """
        Function to set all values of active variables to the values in the large np.array x.
        """
        for i, var in enumerate(self.getActiveVariablesIsWriteable()):
            var.setvalue(x[i])


    def setStatus(self, name, var_status=True):
        """
        Change status of certain variables.
        """
        self.dict_variables[name].status = var_status

def MeritFunctionWrapperScipy(x, s, meritfunction):
    """
    Merit function wrapper for scipy optimize. Notice that x and length of active values must have the same size

    :param x (np.array): active variable values
    :param s (ClassWithOptimizableVariables): system to be optimized
    :param meritfunction (function): meritfunction depending on s

    :return value of the merit function
    """
    s.setActiveVariableValues(x)

    return meritfunction(s)

def optimizeNewton1D(s, meritfunction, dx, iterations=1):
    """
    Optimization function: Newton1D
    """
    pass

def optimizeSciPyInterface(s, meritfunction, **kwargs):
    """
    Optimization function: Scipy.optimize wrapper
    """
    x0 = s.getActiveVariableValues()
    res = minimize(MeritFunctionWrapperScipy, x0, args=(s, meritfunction), method=kwargs["method"])
    print res
    s.setActiveVariableValues(res.x)
    return s

def optimizeSciPyNelderMead(s, meritfunction, **kwargs):
    """
    Optimization function: direct access to Nelder-Mead algorithm in Scipy.
    """
    return optimizeSciPyInterface(s, meritfunction, method="Nelder-Mead")


if __name__ == "__main__":
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


    class ExampleOS(ClassWithOptimizableVariables):
        def __init__(self):
            super(ExampleOS, self).__init__()
            self.addVariable("X", OptimizableVariable(True, "Variable", value=10.0))
            self.addVariable("Y", OptimizableVariable(True, "Variable", value=20.0))
            self.addVariable("Z",
                         OptimizableVariable(True, "Solve",
                                             function=lambda x, y: x**2 + y**2,
                                             args=(self.dict_variables["X"], self.dict_variables["Y"])))


    def testmerit(s):
        return s.dict_variables["X"].evaluate()**2 \
            + s.dict_variables["Y"].evaluate()**2 \
            + s.dict_variables["Z"].evaluate()**2


    def f(p, q):
        return p + q

    p = OptimizableVariable(False, "Variable", value="glass1")
    q = OptimizableVariable(True, "Variable", value="glass2")
    r = OptimizableVariable(False, "Solve", function=f, args=(p, q))
    s = OptimizableVariable(False, "Pickup", function=f, args=(1.0, 6.0))
    print p
    print p.__dict__
    print p.evaluate()
    print q
    print q.__dict__
    print q.evaluate()
    print r
    print r.__dict__
    print r.evaluate()
    print s
    print s.__dict__
    print s.evaluate()
    p.setvalue("glass5") # TODO: assignment operator overloading
    q.setvalue("glass6") # TODO: should behave different for Variable or Solve
    print r.evaluate()
    r.parameters["function"] = lambda x, y: x*3 + y*4
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

    os = ExampleOS()
    print os.dict_variables["X"]
    print os.dict_variables["Y"]

    optimizeSciPyNelderMead(os, testmerit)
    print os.dict_variables["X"]
    print os.dict_variables["Y"]
    print os.dict_variables["Z"].evaluate()

