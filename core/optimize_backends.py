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

from scipy.optimize import minimize
import numpy as np

class Backend(object):
    """
    Base class for the optimization backend. Performs one full optimization run.
    Eats 1D numpy array as starting value. Spits out 1D numpy array as final
    value. Needs to know an optimization function which converts the 1D
    numpy array into the real valued function which is to be optimized.
    """
    
    
    def __init__(self, **kwargs):
        """
        kwargs is everything which is known at initialization time and needed
        by the optimization backend        
        """        
        self.options = kwargs

    def init(self, func):
        """
        Tells backend which function to optimize (usually if coupled to
        optimizer this is a merit function wrapper)
        """
        self.func = func
        
    def run(self, x0):
        """
        Performs optimization. Start value is x0. Has to return xfinal.
        """
        
        raise NotImplementedError()


class ScipyBackend(Backend):
    """
    Uses scipy for optimization.
    """
    
    def run(self, x0):
        res = minimize(self.func, x0, args=(), **self.options)
        return res.x
        
class Newton1DBackend(Backend):
    """
    Uses 1D Newton approach for optimization. Does not need scipy.
    Is only intended as failsafe solution.
    """
    
    
    def run(self, x0):
        
        dx = self.options.get("dx", 1e-6) # set default to 1e-6        
        iters = self.options.get("iterations", 100)
        
        xfinal = x0
        
        for i in range(iters):

            retries = np.ones_like(x0, dtype=bool)        
            retrycount = 0

        
            x0 = xfinal
            while np.all(retries):
                
                retrycount += 1
                
                merit0 = self.func(x0)
                varvalue0 = x0
                varvalue1 = varvalue0 + dx*(1. - 2*np.random.random(np.shape(x0)))
                merit1 = self.func(varvalue1)            
                
                to_be_updated = np.logical_not(np.isclose(merit1 - merit0, 0))
            
                m = (merit1 - merit0) / (varvalue1 - varvalue0)        
                n = merit0 - m * varvalue0
            
                varvalue2 = np.where(to_be_updated, -n/m, varvalue1)
                xfinal = varvalue2
                merit2 = self.func(varvalue2)

                print(n, m, varvalue0, varvalue1, varvalue2)

                
                guard = 0                
                while merit2 > merit0 and np.all(np.abs(varvalue2 - varvalue0)) > dx and guard < 1000:
                    varvalue2 = 0.5*(varvalue2 + varvalue0)
                    xfinal = varvalue2
                    merit2 = self.func(varvalue2)
                    guard += 1
                    
                retries = False
                
        return xfinal
        