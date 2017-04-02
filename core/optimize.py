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
import uuid

class OptimizableVariable(object):
    """
    Class that contains an optimizable variable. Used to get a pointer on a variable.
    The value is not constrained to float. Also other dependent variables are possible to define.
    """
    def __init__(self, variable_status=False, variable_type="Variable", **kwargs):
        """
        Name is gone since it's only needed for reference. Therefore the former listOfOptimizableVariables
        will be a dictionary. type may contain a string e.g.
        "Variable", "Pickup", "External", "...",
        status=True means updating during optimization run
        status=False means no updating during optimization run
        kwargs depend on type
        "Variable" is value=value
        "Pickup" gets function=f, args=tuple of optimizablevariables
        "External" gets function=f, args=tuple of external standard variables (or values)

        :param variable_status (bool): should variable optimized during optimization run?
        :param variable_type (string): which kind of variable do we have? valid choices are: 'variable', 'pickup', 'external'

        :param kwargs (keyword arguments): used to fill up parameters dictionary for variable.

        kwargs:

        type:        kwargs:
        'variable'    value=1.0, or value='hello'
        'pickup'        function=f, args=(other_optvar1, other_optvar2, ...)
        'external'        function=f, args=(other_externalvar/value, ...)

        Notice: f must take exactly as many arguments as args=... long is. The only constraint on f is
        that it has to convert the variables into some final float variable. If this is not the case some
        low-level optimizer might break down.


        """
        self.__var_type = variable_type.lower()
        self.status = variable_status
        self.parameters = kwargs

        # TODO: status variable for "pickup" and "external"
        # TODO: either status variable is important or not
        # TODO: most simple solution -> ignore status variable for "pickup" and "external"
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

    def changetype(self, vtype, **kwargs):
        self.var_type = vtype
        self.parameters = kwargs

    def setvalue(self, value):
        if self.var_type == "variable":
            self.parameters["value"] = value

    def eval_variable(self):
        # if type = variable then give only access to value
        return self.parameters["value"]

    def eval_pickup(self):
        # if type = pickup then pack all arguments into one tuple
        # and put it into the userdefined function
        # evaluate the result
        arguments_for_function_eval = (argfunc.evaluate() for argfunc in self.parameters["args"])
        return self.parameters["function"](*arguments_for_function_eval)

    def eval_external(self):
        # same as for solve except that there are no further OptimizableVariables to be considered
        return self.parameters["function"](*self.parameters["args"])



    def evaluate(self):
        # notice: evaluation code is not limited to floats
        # do not use if-else construct because for many cases dict lookup is faster
        evaldict = {
                    "variable" : self.eval_variable,
                    "pickup" : self.eval_pickup,
                    "external" : self.eval_external
                    }

        return evaldict[self.var_type]()
        
class ClassWithOptimizableVariables(object):
    """
    Implementation of some class with optimizable variables with the help of a dictionary.
    This class is also able to collect the variables and their values from its subclasses per recursion.
    """
    def __init__(self, name = ""):
        """
        Initialize with empty dict.
        """
        self.setName(name)
        self.dict_variables = {}
        self.list_observers = [] 
        # for the optimizable variable class it is useful to have some observer links
        # they get informed if variables change their values

    def setName(self, name):
        if name == "":
            name = str(uuid.uuid4())
        self.__name = name
        
    def getName(self):
        return self.__name
        
    name = property(getName, setName)


    def appendObservers(self, obslist):
        self.list_observers += obslist

    def informObservers(self):
        for obs in self.list_observers:
            obs.informAboutUpdate()

    def addVariable(self, name, var):
        """
        Add some variable into dict.
        """
        self.dict_variables[name] = var


                
    def getAllVariables(self):

        def addOptimizableVariablesToList(var, listOfOptVars = [], idlist=[]):
            """
            Accumulates optimizable variables in var and its linked objects.
            Ignores ring-links and double links.       
            
            @param var: object to evaluate (object)
            @param listOfOptVars: optimizable variables found so far (list of objects)
            @param idlist: ids of objects already evaluated (list of int)
            """ 
            
            if id(var) not in idlist:
                idlist.append(id(var))
    
                if isinstance(var, ClassWithOptimizableVariables):
                    for v in var.__dict__.values():
                        listOfOptVars, idlist = addOptimizableVariablesToList(v, listOfOptVars, idlist) 
                elif isinstance(var, dict):
                    for v in var.values():
                        listOfOptVars, idlist = addOptimizableVariablesToList(v, listOfOptVars, idlist)                        
    
                elif isinstance(var, list) or isinstance(var, tuple):
                    for v in var:
                        listOfOptVars, idlist = addOptimizableVariablesToList(v, listOfOptVars, idlist)
    
                elif isinstance(var, OptimizableVariable):
                    listOfOptVars.append(var)
    
            return listOfOptVars, idlist



        (lst, idlist) = addOptimizableVariablesToList(self)
        return lst


    def getAllVariablesOld(self):
        """
        Conversion of dict into list of Variables. These are only references to the objects in dict.
        Therefore all changes to them directly affects the dictionary in the class.
        """
        # if there are sub classes which are inherited from ClassWith... append their getAllVariables()
        # since these are also containing references to the objects in their respective dicts there should
        # not be any problem that there are copies created

        lst_of_vars = self.dict_variables.values()
        lst_of_attributes_which_are_class_with_opt_vars = \
            filter(lambda x: isinstance(x, ClassWithOptimizableVariables), \
                self.__dict__.values())

        #id vergleich!

        for list_vars in \
            filter(lambda x: isinstance(x, list) or isinstance(x, tuple), \
                self.__dict__.values()):
            for a in list_vars:
                if isinstance(a, ClassWithOptimizableVariables):
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
        Returns a np.array of bools which is True for a variable and False for a pickup or external,
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


class Optimizer(object):
    '''
    Easy optimization interface. All variables are public such that a quick
    attachment of other meritfunctions or other update functions with other
    parameters is possible.
    '''
    def __init__(self, classwithoptvariables, meritfunction, updatefunction, meritparameters={}, updateparameters={}):
        super(Optimizer, self).__init__()
        self.classwithoptvariables = classwithoptvariables
        self.meritfunction = meritfunction
        self.updatefunction = updatefunction
        self.meritparameters = meritparameters
        self.updateparameters = updateparameters
        self.log = "" # for logging

    def MeritFunctionWrapperScipy(self, x): #, s, meritfunction, func):
        """
        Merit function wrapper for scipy optimize. Notice that x and length of active values must have the same size
    
        :param x (np.array): active variable values
        :param s (ClassWithOptimizableVariables): system to be optimized
        :param meritfunction (function): meritfunction depending on s
    
        :return value of the merit function
        """
        self.classwithoptvariables.setActiveVariableValues(x)
        self.updatefunction(self.classwithoptvariables, **self.updateparameters)    
        return self.meritfunction(self.classwithoptvariables, **self.meritparameters)

    def optimizeSciPyInterface(self, **kwargs):
        """
        Optimization function: Scipy.optimize wrapper
        """
        steps = kwargs.get("steps", 0.0)
        if steps > 0:
            opts["maxiter"] = steps
        else:
            opts = {}
        
        x0 = self.classwithoptvariables.getActiveVariableValues()
        print(x0) # TODO: rewrite to log
        res = minimize(self.MeritFunctionWrapperScipy, x0, args=(), method=kwargs["method"], options=opts)
        print res # TODO: rewrite to log
        self.classwithoptvariables.setActiveVariableValues(res.x)
        return self.classwithoptvariables

    def optimizeSciPyNelderMead(self, **kwargs):
        """
        Optimization function: direct access to Nelder-Mead algorithm in Scipy.
        """
        return self.optimizeSciPyInterface(method="Nelder-Mead")

    def run(self, optimizer="optimizeSciPyNelderMead", steps=0, **kwargs):
        '''
        Funtion to perform a certain number of optimization steps.
        '''
        eval("self." + optimizer)(steps=steps, **kwargs)



def optimizeNewton1D(s, meritfunction, dx, iterations=1):
    """
    Optimization function: Newton1D
    """
    pass


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
                         OptimizableVariable(True, "Pickup",
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
    r = OptimizableVariable(False, "Pickup", function=f, args=(p, q))
    s = OptimizableVariable(False, "External", function=f, args=(1.0, 6.0))
    print("print p variable")
    print p
    print p.__dict__
    print p.evaluate()
    print("print q variable")
    print q
    print q.__dict__
    print q.evaluate()
    print("print r variable")
    print r
    print r.__dict__
    print r.evaluate()
    print("print s variable")
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
    
    def optnonupdate(s):
        pass

    optimi = Optimizer(os, testmerit, optnonupdate)
    optimi.optimizeSciPyNelderMead()
    #optimi.run("optimizeSciPyNelderMead", steps=1)
    #optimizeSciPyNelderMead(os, testmerit, function=optnonupdate)
    print os.dict_variables["X"]
    print os.dict_variables["Y"]
    print os.dict_variables["Z"].evaluate()
    
    print("NEW IT FUNCTION1")
    print([v.evaluate() for v in os.getAllVariables()])
    print("NEW IT FUNCTION2")
    print([v.evaluate() for v in cl.getAllVariables()])
    print("NEW IT FUNCTION3")
    print([v.evaluate() for v in cl2.getAllVariables()])


    print(os.name)

    print(os.dict_variables.items())

