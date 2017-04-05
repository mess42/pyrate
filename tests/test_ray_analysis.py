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
from core.ray import RayBundle
from core.ray_analysis import RayBundleAnalysis

def test_centroid():
    
    raybundle = RayBundle(x0 = np.array([[1, 0, 0, 1, 2], [0, 1, 0, 1, 2], [0, 0, 1, 1, 2]]),
                          k0 = np.zeros((3, 5)), Efield0 = np.zeros((3, 5)))
                          
    rayanalysis = RayBundleAnalysis(raybundle)
    centroid = rayanalysis.getCentroidPosition()
    assert np.allclose(centroid, 4./5.)

def test_rmsspotsize():
    
    raybundle = RayBundle(x0 = np.array([[1, 0, 0, 1, 2], [0, 1, 0, 1, 2], [0, 0, 1, 1, 2]]),
                          k0 = np.zeros((3, 5)), Efield0 = np.zeros((3, 5)))
                          
    rayanalysis = RayBundleAnalysis(raybundle)
    rmssize = rayanalysis.getRMSspotSize(np.array([0, 0, 0]))
    
    assert np.isclose(rmssize, math.sqrt(18.0/4.0))

