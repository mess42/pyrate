# -*- coding: utf-8 -*-
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

"""
This code snippet should be deleted when the code is implemented into the
main code base.
"""

import numpy as np
from core.globalconstants import c0, mu0, eps0, standard_wavelength
import math

class Material():
    def __init__(self, eps):
        
        self.eps = eps
    
    def getEpsTensor(self, x, k):
        
        (num_dim, num_points) = np.shape(x)
        
        return np.repeat(self.eps[:, :, np.newaxis], num_points, axis=2)
        
    def calcXi(self, x, n, kpa, wave=standard_wavelength):

        (num_dims, num_pts) = np.shape(kpa)

        k0 = 2.*math.pi/wave
        
        eps = self.getEpsTensor(x, kpa)
        
        """
        {{a1 -> Scalar[eps[i, -i]], 
          a2 -> Scalar[eps[i, i1]*eps[-i1, -i]], 
          a3 -> Scalar[eps[i, i1]*eps[-i1, i2]*eps[-i2, -i]], 
          a4 -> Scalar[kpa[-i]*kpa[i]], 
          a5 -> Scalar[eps[i, i1]*kpa[-i]*kpa[-i1]], 
          a6 -> Scalar[eps[i, i1]*eps[-i1, i2]*kpa[-i]*kpa[-i2]], 
          a7 -> Scalar[eps[i, i1]*nx1[-i]*nx1[-i1]], 
          a8 -> Scalar[eps[i, i1]*kpa[-i1]*nx1[-i]], 
          a9 -> Scalar[eps[i, i1]*kpa[-i]*nx1[-i1]], 
          a10 -> Scalar[kpa[i]*nx1[-i]], 
          a11 -> Scalar[eps[i, i1]*eps[-i1, i2]*kpa[-i2]*nx1[-i]], 
          a12 -> Scalar[eps[i, i1]*eps[-i1, i2]*kpa[-i]*nx1[-i2]], 
          a13 -> Scalar[eps[i, i1]*eps[-i1, i2]*nx1[-i]*nx1[-i2]]}}
          
          poly = a4*a5*omega^2 + (-(a1*a5) + a6)*omega^4 + ((a1^3 - 3*a1*a2 + 2*a3)*omega^6)/6 + 
          ((2*a10*a5 + a4*(a8 + a9))*omega^2 + (a11 + a12 - a1*(a8 + a9))*omega^4)*xi + 
          ((a5 + a4*a7 + 2*a10*(a8 + a9))*omega^2 + (a13 - a1*a7)*omega^4)*xi^2 + 
          (2*a10*a7 + a8 + a9)*omega^2*xi^3 + a7*omega^2*xi^4     
        """

        a1 = np.einsum('ii...', eps)        
        a2 = np.einsum('ij...,ji...', eps, eps)
        a3 = np.einsum('ij...,jk...,ki...', eps, eps, eps)
        a4 = np.einsum('i...,i...', kpa, kpa)
        a5 = np.einsum('ij...,i...,j...', eps, kpa, kpa)
        a6 = np.einsum('ij...,jk...,i...,k...', eps, eps, kpa, kpa)
        a7 = np.einsum('ij...,i...,j...', eps, n, n)        
        a8 = np.einsum('ij...,j...,i...', eps, kpa, n)
        a9 = np.einsum('ij...,i...,j...', eps, kpa, n)
        a11 = np.einsum('ij...,jk...,k...,i...', eps, eps, kpa, n)
        a12 = np.einsum('ij...,jk...,i...,k...', eps, eps, kpa, n)
        a13 = np.einsum('ij...,jk...,i...,k...', eps, eps, n, n)

        #p4 = a7*omegabar**2
        #p3 = (2*a10*a7 + a8 + a9)*omegabar**2
        #p2 = (a5 + a4*a7 + 2*a10*(a8 + a9))*omegabar**2 + (a13 - a1*a7)*omegabar**4
        #p1 = (2*a10*a5 + a4*(a8 + a9))*omegabar**2 + (a11 + a12 - a1*(a8 + a9))*omegabar**4
        #p0 = a4*a5*omegabar**2 + (-a1*a5 + a6)*omegabar**4 + 1./6.*(a1**3 - 3*a1*a2 + 2*a3)*omegabar**6

        p4 = a7
        p3 = a8 + a9
        p2 = a5 + a4*a7 + (a13 - a1*a7)*k0**2
        p1 = a4*p3 + (a11 + a12 - a1*p3)*k0**2
        p0 = a4*a5 + (-a1*a5 + a6)*k0**2 + 1./6.*(a1**3 - 3*a1*a2 + 2*a3)*k0**4
        
        xiarray = np.zeros((4, num_pts), dtype=complex)
        for i in np.arange(num_pts):       
            polycoeffs = [p4[i], p3[i], p2[i], p1[i], p0[i]]
            roots = np.roots(polycoeffs)
            xiarray[:, i] = np.pad(roots, (0, 4 - roots.shape[0]), 'constant', constant_values=(np.nan,))
            
        return xiarray
        
if __name__=="__main__":
    
    m = Material(np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]]))
    
    x = np.zeros((3, 5))
    n = np.zeros((3, 5))
    n[2,:] = 1.
    
    k = np.zeros((3, 5))
    k[2,:] = 2.*math.pi/standard_wavelength
    kinplane = k - np.sum(n*k, axis=0)*n    
    
    xis = m.calcXi(x, n, kinplane)

    print(np.array_str(xis, suppress_small=True))
    print(np.array_str(k[2, :], suppress_small=True))
    
