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
from core.material import ConstantIndexGlass
from core.globalconstants import standard_wavelength

def test_anisotropic_xi_calculation():
    lc = LocalCoordinates("1")
    m = ConstantIndexGlass(lc, 1.0)
    wave = standard_wavelength    
    
    x = np.zeros((3, 5))
    n = np.zeros((3, 5))
    n[2, :] = 1.    
    
    kpa = np.zeros((3, 5))
    phi = np.linspace(0, 360, 5, endpoint=False)*math.pi/180.0
    print(phi)
    k0 = 2.*math.pi/wave

    kpa[0, :] = k0*np.cos(phi)
    kpa[1, :] = k0*np.sin(phi)
    
    # FIXME: determinantenformel mit spuren checken!
    
    print(m.calcXiAnisotropic(x, n, kpa, wave=wave))
    
if __name__=="__main__":
    test_anisotropic_xi_calculation()

