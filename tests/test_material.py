#!/usr/bin/env/python
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
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
"""
import numpy as np
import math
from core.localcoordinates import LocalCoordinates
from core.material import AnisotropicMaterial
from core.globalconstants import standard_wavelength

def test_anisotropic_xi_calculation():
    lc = LocalCoordinates("1")
    
    myeps = np.random.random((3, 3))
    ((epsxx, epsxy, epsxz), (epsyx, epsyy, epsyz), (epszx, epszy, epszz)) = \
        tuple(myeps)
    
    print(myeps)
    print('epsxx: ', epsxx)    
    print('epsxy: ', epsxy)    
    print('epsxz: ', epsxz)    
    print('epsyx: ', epsyx)    
    print('epsyy: ', epsyy)    
    print('epsyz: ', epsyz)    
    print('epszx: ', epszx)    
    print('epszy: ', epszy)    
    print('epszz: ', epszz)    
    
    m = AnisotropicMaterial(lc, myeps)
    
    n = np.zeros((3, 5))
    n[2, :] = 1.    
    
    x = np.zeros((3, 5))
    kpa = np.random.random((3, 5))
    kpa[2, :] = 0.

    (p4v, p3v, p2v, p1v, p0v) = m.calcXiPolynomial(x, kpa, n)


    p4 = epszz
    p3 = (epsxz + epszx)*kpa[0, :] + (epsyz + epszy)*kpa[1, :]
    #p2 = -2*ind**2*(-1 + ind**2)
    #p1 = 0.
    #p0 = ind**2*(-1 + ind**2)**2
    
    print("ana p4", p4)    
    print("code p4", p4v)
   
    print("ana p3", p3)    
    print("code p3", p3v)
    #print("p2", p2)    
    #print("p1", p1)    
    #print("p0", p0)    

  
if __name__=="__main__":
    test_anisotropic_xi_calculation()

