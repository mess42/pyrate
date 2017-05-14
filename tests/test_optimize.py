"""
Pyrate - Optical raytracing based on Python

Copyright (C) 2017 Moritz Esslinger <moritz.esslinger@web.de>
               and Johannes Hartung <j.hartung@gmx.net>
               and     Uwe Lippmann <uwe.lippmann@web.de>
               and    Thomas Heinze <t.heinze@fn.de>

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 2
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
MA  02110-1301, USA.
"""

from core.optimize import ClassWithOptimizableVariables, OptimizableVariable, Optimizer
from core.optimize_backends import ScipyBackend, Newton1DBackend

import numpy as np

#class ExampleSubClass(ClassWithOptimizableVariables):
#    def __init__(self):
#        self.a = 10.
#
#class ExampleSuperClass(ClassWithOptimizableVariables):
#    def __init__(self):
#        super(ExampleSuperClass, self).__init__()
#        self.b = 20.
#        self.c = ClassWithOptimizableVariables()
#        self.c.addVariable("blubberbla", OptimizableVariable("Variable", value=5.0))
#        self.addVariable("blubberdieblub", OptimizableVariable("Variable", value=10.0))




def test_variables_pickups_externals():
    """
    Check whether pickups are also working for strings
    """
    def f(p, q):
        return p + q

    p = OptimizableVariable("Variable", value="glass1")
    q = OptimizableVariable("Variable", value="glass2")
    r = OptimizableVariable("Pickup", function=f, args=(p, q))
    s = OptimizableVariable("External", function=f, args=(1.0, 6.0))

    assert p() == "glass1"
    assert q() == "glass2"
    assert r() == "glass1glass2"
    assert s() == 7.0

    p.setvalue("glass5") 
    q.setvalue("glass6") 
    
    assert r() == "glass5glass6"
    
    r.parameters["function"] = lambda x, y: x*3 + y*4

    assert r() == "glass5glass5glass5glass6glass6glass6glass6"


def test_optimization():

    class ExampleOS(ClassWithOptimizableVariables):
        def __init__(self):
            super(ExampleOS, self).__init__()
            self.addVariable("X", OptimizableVariable("Variable", value=3.0))
            self.addVariable("Y", OptimizableVariable("Variable", value=20.0))
            self.addVariable("Z",
                         OptimizableVariable("Pickup",
                                             function=lambda x, y: x**2 + y**2,
                                             args=(self("X"), self("Y"))))
    
    
    def testmerit(s):
        # let x, y being on a circle of radius 5
        return (s("X")()**2 + s("Y")()**2 - 5.**2)**2
        
            

    def testmerit2(s):
        # let x**2 + y**2 + z**2 be minimized with the pickup constraint z = x**2 + y**2
        return s("X")()**2 + s("Y")()**2 + s("Z")()**2
        

    os = ExampleOS()

    optimi = Optimizer(os, testmerit, backend=ScipyBackend(method='Nelder-Mead', options={'maxiter':1000, 'disp':False}))
    optimi.run()
    assert np.isclose(os("X")()**2 + os("Y")()**2, 25.0)
    optimi.meritfunction = testmerit2
    optimi.run()
    assert np.isclose(os("X")()**2 + os("Y")()**2, os("Z")())
    optimi.backend = Newton1DBackend(dx=1e-6, iterations=500)
    optimi.meritfunction = testmerit
    optimi.run()
    assert np.isclose(os("X")()**2 + os("Y")()**2, 25.0)
    optimi.meritfunction = testmerit2
    optimi.run()
    assert np.isclose(os("X")()**2 + os("Y")()**2, os("Z")())

    
