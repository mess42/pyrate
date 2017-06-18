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

import math
import numpy as np
import core.helpers_math
from core.globalconstants import canonical_ex, canonical_ey, standard_wavelength
from core.material import AnisotropicMaterial
from core.localcoordinates import LocalCoordinates
from core.helpers import build_pilotbundle2
from core.surface import Surface        

if __name__=="__main__":
    
    num_pts = 1    
    #rnd_vecs = np.random.randn(3, num_pts)

    rnd_vecs = np.zeros((3, num_pts))
    rnd_vecs[2, :] = 1

    rnd_units = rnd_vecs/np.linalg.norm(rnd_vecs, axis=0)

    #rnd_data1 = np.random.random((3, 3)) #np.eye(3)
    #rnd_data2 = np.random.random((3, 3))#np.zeros((3, 3))#

    rnd_data1 = np.eye(3)
    rnd_data2 = np.zeros((3, 3))

    lc = LocalCoordinates("1")
    myeps = rnd_data1 + complex(0, 1)*rnd_data2

    mat = AnisotropicMaterial(lc, myeps)
    surfobj = Surface(lc)

    pilotbundles = build_pilotbundle2(surfobj, mat, (0.1, 0.1), (0.05, 0.05), kunitvector=None, wave=standard_wavelength)

    fig = plt.figure(1)
    ax = fig.add_subplot(111)
    ax.axis('equal')
    
    pilotbundles[0].draw2d(ax, color="blue")