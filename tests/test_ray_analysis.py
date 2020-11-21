#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Pyrate - Optical raytracing based on Python

Copyright (C) 2014-2020
               by     Moritz Esslinger moritz.esslinger@web.de
               and    Johannes Hartung j.hartung@gmx.net
               and    Uwe Lippmann  uwe.lippmann@web.de
               and    Thomas Heinze t.heinze@uni-jena.de
               and    others

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
from pyrateoptics.raytracer.ray import RayBundle
from pyrateoptics.raytracer.analysis.ray_analysis import RayBundleAnalysis


def test_centroid():
    """
    TODO
    """
    raybundle = RayBundle(x0=np.array([[1, 0, 0, 1, 2],
                                       [0, 1, 0, 1, 2],
                                       [0, 0, 1, 1, 2]]),
                          k0=np.zeros((3, 5)), Efield0=np.zeros((3, 5)))
    rayanalysis = RayBundleAnalysis(raybundle)
    centroid = rayanalysis.get_centroid_position()
    assert np.allclose(centroid, 4./5.)


def test_rmsspotsize():
    """
    TODO
    """
    raybundle = RayBundle(x0=np.array([[1, 0, 0, 1, 2],
                                       [0, 1, 0, 1, 2],
                                       [0, 0, 1, 1, 2]]),
                          k0=np.zeros((3, 5)), Efield0=np.zeros((3, 5)))
    rayanalysis = RayBundleAnalysis(raybundle)
    rmssize = rayanalysis.get_rms_spot_size(np.array([0, 0, 0]))
    assert np.isclose(rmssize, math.sqrt(18.0/4.0))


def test_arc_length():
    """
    TODO
    """
    k0 = np.zeros((3, 2))
    E0 = np.zeros((3, 2))
    raybundle = RayBundle(x0=np.array([[0, 0], [0, 0], [0, 0]]),
                          k0=k0, Efield0=E0)
    x1 = np.array([[1, 0], [0, 0], [0, 0]])
    x2 = np.array([[1, 1], [1, 1], [0, 0]])
    x3 = np.array([[0, 2], [1, 2], [0, 0]])
    x4 = np.array([[0, 3], [0, 3], [0, 0]])
    valid = np.ones_like([1, 1])
    raybundle.append(x1, k0, E0, valid)
    raybundle.append(x2, k0, E0, valid)
    raybundle.append(x3, k0, E0, valid)
    raybundle.append(x4, k0, E0, valid)
    arclen = RayBundleAnalysis(raybundle).get_arc_length()
    assert np.allclose(arclen, np.array([4., 3*np.sqrt(2)]))

def test_direction_centroid():
    """
    TODO
    """
    k0 = np.zeros((3, 5))
    k0[2, :] = 1
    E0 = np.zeros((3, 5))
    E0[1, :] = 1.
    raybundle = RayBundle(x0=np.array([[1, 0, 0, 1, 2],
                                       [0, 1, 0, 1, 2],
                                       [0, 0, 1, 1, 2]]),
                          k0=k0, Efield0=E0)
    rayanalysis = RayBundleAnalysis(raybundle)
    centroiddir = rayanalysis.get_centroid_direction()
    assert np.allclose(centroiddir, np.array([0, 0, 1]))


def test_rms_angularsize():
    """
    TODO
    """
    k0 = np.zeros((3, 5))
    k0[2, :] = 1
    E0 = np.zeros((3, 5))
    E0[1, :] = 1.
    raybundle = RayBundle(x0=np.array([[1, 0, 0, 1, 2],
                                       [0, 1, 0, 1, 2],
                                       [0, 0, 1, 1, 2]]),
                          k0=k0, Efield0=E0)
    rayanalysis = RayBundleAnalysis(raybundle)
    angularsize = rayanalysis.get_rms_angluar_size(
        np.array([math.sin(1.*math.pi/180.0), 0, math.cos(1.*math.pi/180.0)]))
    assert np.isclose(angularsize, (1.*math.pi/180.0))
