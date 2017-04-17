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

from core.optimize import ClassWithOptimizableVariables, OptimizableVariable

class CoreObject(ClassWithOptimizableVariables):
    def __init__(self):
        super(CoreObject, self).__init__()
        self.addVariable("X", OptimizableVariable("Variable", value=3.0))
        self.addVariable("Y", OptimizableVariable("Variable", value=20.0))
        self.addVariable("Z",
                     OptimizableVariable("Pickup",
                                         function=lambda x, y: x**2 + y**2,
                                         args=(self.dict_variables["X"], self.dict_variables["Y"])))

    def optimize(self, optimizer):
        optimizer.visit(self)
        # optimizer.visit(self.dependenctObject)   
        # the model knows best
        # about its dependencies

class Optimizer:

    def __init__(self):    
        pass

    def optimizeCoreObject(self, coreobj, meritfunc, updatefunc):
        # For every core object there is a certain user level merit function
        # necessary. Also the updatefunc is necessary, to e.g. update subsystems
        # during optimization run
        pass
       
    def visit(self, obj):
        # how to call optimizeCoreObject?
        pass


class ScipyOptimizer(Optimizer):
        
    def optimizeCoreObject(self, coreobj, meritfunc, updatefunc):
        # do optimizations, where SciPy is used only here
        pass


if __name__ == "__main__":
    co = CoreObject()

    def updatefunc(coreobj):
        pass
    
    def meritfunc(coreobj):
        return (coreobj.dict_variables["X"].evaluate()**2 \
            + coreobj.dict_variables["Y"].evaluate()**2 - 5.**2)**2#\

    print(meritfunc(co))    
    
    so = ScipyOptimizer()
    
    co.optimize(so)