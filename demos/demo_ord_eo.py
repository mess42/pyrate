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

"""
Smoke demo for picture checking. This means the result of this file
is compared with a fixed picture and may not change. That's part of the test.
If anything breaks, the picture  usually changes.
"""


import numpy as np
import matplotlib.pyplot as plt

import matplotlib
from distutils.version import StrictVersion


from core import raster
from core import material
from core import surfShape
from core.optical_element import OpticalElement
from core.surface import Surface
from core.optical_system import OpticalSystem
from core.ray import RayBundle

from core.aperture import CircularAperture
from core.localcoordinates import LocalCoordinates

from core.globalconstants import canonical_ey

import math

wavelength = 0.5876e-3

# definition of optical system
s = OpticalSystem() 

lc0 = s.addLocalCoordinateSystem(LocalCoordinates(name="stop", decz=1.0), refname=s.rootcoordinatesystem.name)
lc1 = s.addLocalCoordinateSystem(LocalCoordinates(name="surf1", decz=10.0), refname=lc0.name) # objectDist
lc2 = s.addLocalCoordinateSystem(LocalCoordinates(name="surf2", decz=5.0), refname=lc1.name)
lc3 = s.addLocalCoordinateSystem(LocalCoordinates(name="image", decz=10.0), refname=lc2.name)


stopsurf = Surface(lc0)
frontsurf = Surface(lc1, shape=surfShape.Conic(lc1, curv=0), apert=CircularAperture(lc1, 10.0))
rearsurf = Surface(lc2, shape=surfShape.Conic(lc2, curv=0), apert=CircularAperture(lc3, 10.0))
image = Surface(lc3)


elem = OpticalElement(lc0, label="crystalelem")

no = 1.5
neo = 1.8

myeps = np.array([[no, 0, 0], [0, no, 0], [0, 0, neo]]) 

crystal = material.AnisotropicMaterial(lc1, myeps)


elem.addMaterial("crystal", crystal)

elem.addSurface("stop", stopsurf, (None, None))
elem.addSurface("front", frontsurf, (None, "crystal"))
elem.addSurface("rear", rearsurf, ("crystal", None))
elem.addSurface("image", image, (None, None))

s.addElement("crystalelem", elem)

rstobj = raster.MeridionalFan()
(px, py) = rstobj.getGrid(10)

rpup = 8.0

phik = 20.*math.pi/180.0

o = np.vstack((np.zeros_like(px), np.zeros_like(px), -5*np.ones_like(px)))
k = np.zeros_like(o)
k0 = 2.*math.pi/wavelength
k[1, :] = k0*np.sin(phik*py)
k[2, :] = k0*np.cos(phik*py)
#o = np.vstack((rpup*px, rpup*py, -5.*np.ones_like(px)))
#k = np.zeros_like(o)
#k[1,:] = k0*math.sin(phik)
#k[2,:] = k0*math.cos(phik)

ey = np.zeros_like(o)
ey[1,:] =  1.

E0 = np.cross(k, ey, axisa=0, axisb=0).T

sysseq = [("crystalelem", [("stop", True, True), ("front", True, True), ("rear", True, True), ("image", True, True)])]

phi = 5.*math.pi/180.0

initialbundle = RayBundle(x0=o, k0=k, Efield0=E0, wave=wavelength)
r2 = s.seqtrace(initialbundle, sysseq, splitup=False)

fig = plt.figure(1)
ax = fig.add_subplot(111)

ax.axis('equal')
if StrictVersion(matplotlib.__version__) < StrictVersion('2.0.0'):
    ax.set_axis_bgcolor('white')
else:
    ax.set_facecolor('white')

phi = 0.#math.pi/4
pn = np.array([math.cos(phi), 0, math.sin(phi)]) # canonical_ex
up = canonical_ey

r2[0].draw2d(ax, color="blue", plane_normal=pn, up=up) 
#r2[1].draw2d(ax, color="green", plane_normal=pn, up=up) 
for e in s.elements.itervalues():
    for surfs in e.surfaces.itervalues():
        surfs.draw2d(ax, color="grey", vertices=50, plane_normal=pn, up=up) # try for phi=0.
        #surfs.draw2d(ax, color="grey", inyzplane=False, vertices=50, plane_normal=pn, up=up) # try for phi=pi/4


plt.show()


