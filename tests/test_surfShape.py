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
import math
from core.surfShape import Conic, Asphere
from core.localcoordinates import LocalCoordinates

def test_conic_sag():
    
    lc = LocalCoordinates(name="root")    
    
    radius = 10.0
    cc = -1.5
    curv = 1./radius
        
    if cc > -1:
        maxradius = math.sqrt(1./((1+cc)*curv**2))    
    else:
        maxradius = radius
        
    pts = (2*np.random.random((2, 10))- 1.)*maxradius    
    x = pts[0]
    y = pts[1]
    
    shape = Conic(lc, curv=curv, cc=cc)    
    z = shape.getSag(x, y)    
    # comparison with explicitly entered formula
    assert np.allclose(z, curv*(x**2 + y**2)/(1. + np.sqrt(1. - (1.+cc)*curv**2*(x**2 + y**2))))

def test_asphere_sag():
    
    radius = 10.0
    cc = 0
    a2 = 1e-3
    a4 = -1e-6
    a6 = 1e-8

    lc = LocalCoordinates(name="root")    
    
    radius = 10.0
    cc = -1.5
    curv = 1./radius
        
    if cc > -1:
        maxradius = math.sqrt(1./((1+cc)*curv**2))    
    else:
        maxradius = radius
        
    pts = (2*np.random.random((2, 10))- 1.)*maxradius    
    x = pts[0]
    y = pts[1]

    
    shape = Asphere(lc, curv=curv, cc=cc, acoeffs=[a2, a4, a6])    
    z = shape.getSag(x, y)
    # comparison with explicitly entered formula
    comparison = curv*(x**2 + y**2)/(1. + np.sqrt(1. - (1.+cc)*curv**2*(x**2 + y**2))) \
                + a2*(x**2 + y**2) + a4*(x**2 + y**2)**2 + a6*(x**2 + y**2)**3
    assert np.allclose(z, comparison)


#if __name__=="__main__":
#    test_conic()
#    test_asphere()