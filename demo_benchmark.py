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
import matplotlib.pyplot as plt
import time
import math

from core import pupil
from core import raster
from core import material
from core import surfShape
from core.optical_system import OpticalSystem
from core.ray import RayPath, RayBundle

from core import plots
from core.aperture import CircularAperture
from core.localcoordinates import LocalCoordinates

from core.globalconstants import standard_wavelength

from core.helpers import build_simple_optical_system

wavelength = standard_wavelength

# definition of optical system

mat_dict = {"glass":1.7, "glass2":1.5}


s = build_simple_optical_system(
        [(-5.922, 0, 2.0, None, ""),
         (-3.160, 0, 3.0, "glass", ""),
         (15.884, 0, 5.0, None, ""),
        (-12.756, 0, 3.0, "glass", ""),
        (0, 0, 3.0, None, ""),
        (3.125, 0, 2.0, "glass2", ""),
        (1.479, 0, 3.0, None, ""),
        (0, 0, 19.0, None, "")
         ], mat_dict)


# benchmark

# pilot bundle

pts = np.array([[0,0, 0], [0.1, 0.2, 0.3], [0, 0, 0]])
dirs = np.array([[0,0, 0], [0, 0, 0], [1, 1, 1]])

pilotbundle = RayBundle(pts, dirs, s.surfaces[0].material, np.array([0, 1, 2]), wave=wavelength, pol=[])
pilotpath = RayPath(pilotbundle, s)

#print([blub.o for blub in pilotpath.raybundles])

# definition of rays
nray = 1E5 # number of rays
aimy = aim.aimFiniteByMakingASurfaceTheStop(s, pupilType=pupil.ObjectSpaceNA, #.StopDiameter,
                                            pupilSizeParameter=0.2,#3.0,
                                            fieldType= field.ObjectHeight,
                                            rasterType= raster.RectGrid,
                                            nray=nray, wavelength=wavelength, stopPosition=5)
initialBundle = aimy.getInitialRayBundle(s, fieldXY=np.array([0, 0]), wavelength=wavelength)
nray = len(initialBundle.o[0, :])

t0 = time.clock()
r = RayPath(initialBundle, s)
print "benchmark : ", time.clock() - t0, "s for tracing ", nray, " rays through ", len(s.surfaces) - 1, " surfaces."
print "             That is ", int(round(nray * (len(s.surfaces) - 1) / (time.clock() - t0))), "ray-surface-operations per second"

# plot
aimy.setPupilRaster(rasterType= raster.RectGrid, nray=100)

initialBundle2 = aimy.getInitialRayBundle(s, fieldXY=np.array([0, 0]), wavelength=wavelength)
r2 = RayPath(initialBundle2, s)

initialBundle3 = aimy.getInitialRayBundle(s, fieldXY=np.array([0, 0.1]), wavelength=wavelength)
r3 = RayPath(initialBundle3, s)

fig = plt.figure(1)
ax = fig.add_subplot(111)

ax.axis('equal')
ax.set_axis_bgcolor('black')

#plots.drawLayout2d(ax, s, [pilotpath])
plots.drawLayout2d(ax, s, [r2,r3])

print(s.getParaxialPupil( stopPosition=5, ray=initialBundle))

plt.show()


