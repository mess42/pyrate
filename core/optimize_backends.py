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
            
            #print("iteration %d merit value: %f" % (i, self.func(x0)))            
            
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

                guard = 0                
                while merit2 > merit0 and np.all(np.abs(varvalue2 - varvalue0)) > dx and guard < 1000:
                    varvalue2 = 0.5*(varvalue2 + varvalue0)
                    xfinal = varvalue2
                    merit2 = self.func(varvalue2)
                    guard += 1
                    
                retries = False
                
        return xfinal
        
class ParticleSwarmBackend(Backend):
        
    def run(self, x0):

        class particle:
            
            def __init__(self, x0, v0):
                self.x = x0
                self.v = v0
                self.pb = x0
        
        initcube = self.options.get("cube", np.vstack((-np.ones(np.shape(x0)), np.ones(np.shape(x0)))))
        cubedelta = initcube[1] - initcube[0]
        num_particles = self.options.get("num_particles", 10)
        max_velocities = self.options.get("max_velocities", 0.1*np.ones(np.shape(x0)))
        tol = self.options.get("tol", 1e-3)
        max_iters = self.options.get("iterations", 100)
        
        particle_list = [
            particle(x0 + initcube[0] + np.random.random(np.shape(x0))*cubedelta,
                     max_velocities*(1. - 2.*np.random.random(np.shape(x0)))) for i in range(num_particles)]        

        termination = False
        
        def bestglobalpos():
            result = np.copy(particle_list[0].x)
            for p in particle_list:
                if self.func(p.x) < self.func(result):
                    result = np.copy(p.x)
            return result
            
       
        c1 = self.options.get("c1", 2.0)
        c2 = self.options.get("c2", 2.0)        
        
        iters = 0        
        
        result = np.zeros(np.shape(x0))
        
        
        while not termination and iters < max_iters:
            iters += 1
            pg = bestglobalpos()

            phi = c1 + c2
            chi =  2./np.abs(2. - phi - np.sqrt(phi**2 - 4.*phi))   


            particle_com = np.zeros(np.shape(x0))

            
            for p in particle_list:
                if self.func(p.x) < self.func(p.pb):
                    p.pb = p.x
                    
                r1 = np.random.random()
                r2 = np.random.random()
                
                
                p.v = chi*(p.v + c1*r1*(p.pb - p.x) + c2*r2*(pg - p.x))
                p.x = p.x + p.v
           
                particle_com += p.x/num_particles

            particle_rms = 0
            for p in particle_list:
                particle_rms += np.sum((p.x - particle_com)**2)/num_particles
            
            particle_rms = np.sqrt(particle_rms)
            
            termination = particle_rms < tol
            
            result = pg
        
        return result
        
            

        
if __name__=="__main__":

    def fun(x):
        return ((x[0] - 1)**2 + (x[1] - 2)**2 - 25)**2
    
    p = ParticleSwarmBackend(cube=np.array([[-10, -10], [10, 10]]), c1=2.1, c2=2.1)
    p.func = fun
    
    
    xfinal = p.run(np.array([20, 20]))
    
    print(xfinal)