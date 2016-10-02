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


from core import pupil
from core import field
from core import raster
from core import material
from core import aim
from core import surfShape
from core.optical_system import OpticalSystem, Surface
from core.ray import RayPath, RayBundle

from core import plots
from core.aperture import CircularAperture
from core.coordinates import LocalCoordinates

import math

# definition of optical system
s = OpticalSystem() # objectDistance = 2.0

lc1 = s.addLocalCoordinateSystem(LocalCoordinates(name="surf1", decz=2.0)) # objectDist
lc2 = s.addLocalCoordinateSystem(LocalCoordinates(name="surf2", decz=3.0))
lc3 = s.addLocalCoordinateSystem(LocalCoordinates(name="surf3", decz=5.0, tiltx=0.0*math.pi/180.0))
lc4 = s.addLocalCoordinateSystem(LocalCoordinates(name="surf4", decz=3.0))
lc5 = s.addLocalCoordinateSystem(LocalCoordinates(name="surf5", decz=3.0))
lc6 = s.addLocalCoordinateSystem(LocalCoordinates(name="surf6", decz=2.0))
lc7 = s.addLocalCoordinateSystem(LocalCoordinates(name="surf7", decz=3.0))
lc8 = s.addLocalCoordinateSystem(LocalCoordinates(name="image", decz=19.0))


s.insertSurface(1, Surface(lc1, surfShape.Conic(curv=1/-5.922), # thickness=3.0,
                           material=material.ConstantIndexGlass(1.7), 
                            aperture=CircularAperture(0.55)))

s.insertSurface(2, Surface(lc2, surfShape.Conic(curv=1/-3.160), # thickness=5.0, 
                           aperture=CircularAperture(1.0)))

s.insertSurface(3, Surface(lc3, surfShape.Conic(curv=1/15.884), #thickness=3.0,
                           material=material.ConstantIndexGlass(1.7), 
                            aperture=CircularAperture(1.3)))

s.insertSurface(4, Surface(lc4, surfShape.Conic(curv=1/-12.756), #thickness=3.0,
                           aperture=CircularAperture(1.3)))

#s.insertSurface(5, Surface(surfShape.Decenter(dx = 0., dy = 1.), material=material.Tilt(angle=20.*np.pi/180.0, axis='X')))

s.insertSurface(5, Surface(lc5, surfShape.Conic(), #thickness=2.0, 
                           aperture=CircularAperture(1.01))) # Stop Surface

s.insertSurface(6, Surface(lc6, surfShape.Conic(curv=1/3.125), #thickness=3.0,
                           material=material.ConstantIndexGlass(1.5), 
                            aperture=CircularAperture(1.0)))

s.insertSurface(7, Surface(lc7, surfShape.Conic(curv=1/1.479), #thickness=19.0,
                           aperture=CircularAperture(1.0)))


s.insertSurface(8, Surface(lc8)) # image


# benchmark

# pilot bundle

pts = np.array([[0,0, 0], [0.1, 0.2, 0.3], [0, 0, 0]])
dirs = np.array([[0,0, 0], [0, 0, 0], [1, 1, 1]])

pilotbundle = RayBundle(pts, dirs, s.surfaces[0].material, np.array([0, 1, 2]), wave=0.55, pol=[])
pilotpath = RayPath(pilotbundle, s)

print([blub.o for blub in pilotpath.raybundles])

# definition of rays
nray = 1E5 # number of rays
aimy = aim.aimFiniteByMakingASurfaceTheStop(s, pupilType=pupil.ObjectSpaceNA, #.StopDiameter,
                                            pupilSizeParameter=0.2,#3.0,
                                            fieldType= field.ObjectHeight,
                                            rasterType= raster.RectGrid,
                                            nray=nray, wavelength=0.55, stopPosition=5)
initialBundle = aimy.getInitialRayBundle(s, fieldXY=np.array([0, 0]), wavelength=.55)
nray = len(initialBundle.o[0, :])

t0 = time.clock()
r = RayPath(initialBundle, s)
print "benchmark : ", time.clock() - t0, "s for tracing ", nray, " rays through ", len(s.surfaces) - 1, " surfaces."
print "             That is ", int(round(nray * (len(s.surfaces) - 1) / (time.clock() - t0))), "ray-surface-operations per second"

# plot
aimy.setPupilRaster(rasterType= raster.RectGrid, nray=100)

initialBundle2 = aimy.getInitialRayBundle(s, fieldXY=np.array([0, 0]), wavelength=.55)
r2 = RayPath(initialBundle2, s)

initialBundle3 = aimy.getInitialRayBundle(s, fieldXY=np.array([0, 0.1]), wavelength=.55)
r3 = RayPath(initialBundle3, s)

fig = plt.figure(1)
ax = fig.add_subplot(111)

ax.axis('equal')
ax.set_axis_bgcolor('black')

#plots.drawLayout2d(ax, s, [pilotpath])
plots.drawLayout2d(ax, s, [r2])

plt.show()

print(s.globalcoordinatesystem.pprint())
