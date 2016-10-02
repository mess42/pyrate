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
import math

from core import pupil
from core import field
from core import raster
from core import material
from core import aim
from core import surfShape
from core import merit
from core import optimize
from core.optical_system import OpticalSystem, Surface
from core.ray import RayPath

from core import plots
from core.aperture import CircularAperture
from core.coordinates import LocalCoordinates

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

s.insertSurface(7, Surface(lc7, surfShape.Conic(curv=0.1*1/1.479), #thickness=19.0,
                           aperture=CircularAperture(1.0)))


s.insertSurface(8, Surface(lc8)) # image




# plot initial system
aimy = aim.aimFiniteByMakingASurfaceTheStop(s, pupilType=pupil.ObjectSpaceNA, #.StopDiameter,
                                            pupilSizeParameter=0.2,#3.0,
                                            fieldType= field.ObjectHeight,
                                            rasterType= raster.RectGrid,
                                            nray=20, wavelength=0.55, stopPosition=5)

initialBundle2 = aimy.getInitialRayBundle(s, fieldXY=np.array([0, 0]), wavelength=.55)

r2 = RayPath(initialBundle2, s)

initialBundle3 = aimy.getInitialRayBundle(s, fieldXY=np.array([0, 0.1]), wavelength=.55)
r3 = RayPath(initialBundle3, s)

fig = plt.figure(1)
ax = fig.add_subplot(311)
ax2 = fig.add_subplot(312)

ax.axis('equal')
ax.set_axis_bgcolor('black')
ax2.axis('equal')
ax2.set_axis_bgcolor('black')

plots.drawLayout2d(ax, s, [r2, r3])


# optimize
print "Initial   merit function: ", merit.mySimpleDumpRMSSpotSizeMeritFunction(s)

# make surface curvatures variable
#s.surfaces[2].shape.setStatus("curvature", True)
#s.surfaces[3].shape.setStatus("curvature", True)
#s.surfaces[4].shape.setStatus("curvature", True)
#s.surfaces[5].shape.setStatus("curvature", True)
#s.surfaces[7].shape.setStatus("curvature", True)


s.surfaces[2].shape.curvature.status=True
s.surfaces[3].shape.curvature.status=(True)
s.surfaces[4].shape.curvature.status=(True)
s.surfaces[5].shape.curvature.status=(True)
s.surfaces[7].shape.curvature.status=(True)


print "aimy,stopDiameter before: ", aimy.stopDiameter

s = optimize.optimizeSciPyInterface(s, merit.mySimpleDumpRMSSpotSizeMeritFunction, method='nelder-mead', options={'xtol': 1e-8, 'disp': True})

print "aimy,stopDiameter after: ", aimy.stopDiameter

print "Optimized merit function: ", merit.mySimpleDumpRMSSpotSizeMeritFunction(s)

aimy.setPupilRaster(rasterType= raster.RectGrid, nray=100)
initialBundle2 = aimy.getInitialRayBundle(s, fieldXY=np.array([0, 0]), wavelength=.55)
r2 = RayPath(initialBundle2, s)
initialBundle3 = aimy.getInitialRayBundle(s, fieldXY=np.array([0, 0.1]), wavelength=.55)
r3 = RayPath(initialBundle3, s)

fig15 = plt.figure(15)
ax3 = fig15.add_subplot(111)

plots.drawLayout2d(ax2, s, [r2, r3])
plots.drawSpotDiagram(ax3, s, r3, -1)

plt.show()



