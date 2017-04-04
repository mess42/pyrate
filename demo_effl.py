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
from core.aperture import CircularAperture, BaseAperture
from core.coordinates import LocalCoordinates

import math

# formula for effective focal length

def effl(r, alpha, phi):
    sp = math.sin(phi)
    cp = math.cos(phi)
    return r*math.sqrt((1. - alpha**2*sp**2)**3/(- 1 + alpha**2*sp**2 + alpha*cp*math.sqrt(1. - alpha**2*sp**2))**2)

def effl_pt(r, alpha, phi):
    sp = math.sin(phi)
    cp = math.cos(phi)

    def Power(x,y):
        return x**y

    def Sqrt(x):
        return math.sqrt(x)

    return np.array([
        -((Power(1 - Power(alpha,2) + Power(alpha,2)*Power(cp,2),2)*r)/(-1 + Power(alpha,2)*Power(sp,2) +
        alpha*cp*Sqrt(1 - Power(alpha,2)*Power(sp,2)))),
        -((alpha*Power(1 - Power(alpha,2) + Power(alpha,2)*Power(cp,2),1.5)*r*sp)/(-1 + Power(alpha,2)*Power(sp,2) +
        alpha*cp*Sqrt(1 - Power(alpha,2)*Power(sp,2))))
                    ]
    )

# definition of optical system

s = OpticalSystem()


lc1 = s.addLocalCoordinateSystem(LocalCoordinates(name="surf1", decz=0.0)) # objectDist
lc2 = s.addLocalCoordinateSystem(LocalCoordinates(name="surf2", decz=10.0))

s.insertSurface(1, Surface(lc1, surfShape.Conic(curv=1./3.),
                           material=material.ConstantIndexGlass(1.7), aperture=CircularAperture(3.0)))
s.insertSurface(2, Surface(lc2))
# pilot bundle

fig = plt.figure(1)
ax = fig.add_subplot(111)

ax.axis('equal')
ax.set_axis_bgcolor('black')


phirange = np.linspace(0.0, 1.0, 5)

phirange_anal = list(np.linspace(0.0, 1.5, 100))

effls = np.zeros((np.shape(phirange_anal)[0], 2), dtype=float)

for (ind, phiangle) in enumerate(phirange):

    pts = np.array([[0, 0], [0.0, 0.1], [0, 0]])
    dirs = np.array([[0,0], np.sin([phiangle, phiangle]), np.cos([phiangle, phiangle])])

    pilotbundle = RayBundle(pts, dirs, s.surfaces[0].material, np.array([0, 1, 2]), wave=standard_wavelength, pol=[])
    pilotpath = RayPath(pilotbundle, s)



    plots.drawLayout2d(ax, s, [pilotpath])

for (ind, blub) in enumerate(phirange_anal):
    effls[ind] = effl_pt(3.0, 1/1.7, blub)

ax.plot(effls[:, 0], effls[:, 1], color='r')

print(effl_pt(3.0, 1/1.7, 1.0))

# definition of rays
#nray = 1E5 # number of rays
#aimy = aim.aimFiniteByMakingASurfaceTheStop(s, pupilType=pupil.ObjectSpaceNA, #.StopDiameter,
#                                            pupilSizeParameter=0.2,#3.0,
#                                            fieldType= field.ObjectHeight,
#                                            rasterType= raster.RectGrid,
#                                            nray=nray, wavelength=0.55, stopPosition=1)
#initialBundle = aimy.getInitialRayBundle(s, fieldXY=np.array([0, 0]), wavelength=.55)
#nray = len(initialBundle.o[0, :])

#t0 = time.clock()
#r = RayPath(initialBundle, s)
#print "benchmark : ", time.clock() - t0, "s for tracing ", nray, " rays through ", len(s.surfaces) - 1, " surfaces."
#print "             That is ", int(round(nray * (len(s.surfaces) - 1) / (time.clock() - t0))), "ray-surface-operations per second"

# plot
#aimy.setPupilRaster(rasterType= raster.RectGrid, nray=20)

#initialBundle2 = aimy.getInitialRayBundle(s, fieldXY=np.array([0, 0]), wavelength=.55)
#r2 = RayPath(initialBundle2, s)

#initialBundle3 = aimy.getInitialRayBundle(s, fieldXY=np.array([0, 0.1]), wavelength=.55)
#r3 = RayPath(initialBundle3, s)



plt.show()



