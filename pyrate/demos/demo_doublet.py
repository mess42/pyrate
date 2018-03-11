#!/usr/bin/env/python
"""
Pyrate - Optical raytracing based on Python

Copyright (C) 2014-2018
               by     Moritz Esslinger moritz.esslinger@web.de
               and    Johannes Hartung j.hartung@gmx.net
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

import matplotlib
from distutils.version import StrictVersion


from core import raster
from core.material_isotropic import ConstantIndexGlass
from core import surfShape
from core.optical_element import OpticalElement
from core.surface import Surface
from core.optical_system import OpticalSystem
from core.ray import RayBundle

from core.aperture import CircularAperture
from core.localcoordinates import LocalCoordinates

from core.globalconstants import canonical_ey

import math
import logging
logging.basicConfig(level=logging.DEBUG)

wavelength = 0.5876e-3

# definition of optical system
s = OpticalSystem() 

lc0 = s.addLocalCoordinateSystem(LocalCoordinates(name="stop", decz=0.0), refname=s.rootcoordinatesystem.name)
lc1 = s.addLocalCoordinateSystem(LocalCoordinates(name="surf1", decz=-1.048), refname=lc0.name) # objectDist
lc2 = s.addLocalCoordinateSystem(LocalCoordinates(name="surf2", decz=4.0), refname=lc1.name)
lc3 = s.addLocalCoordinateSystem(LocalCoordinates(name="surf3", decz=2.5), refname=lc2.name)
lc4 = s.addLocalCoordinateSystem(LocalCoordinates(name="image", decz=97.2), refname=lc3.name)


stopsurf = Surface(lc0)
frontsurf = Surface(lc1, shape=surfShape.Conic(lc1, curv=1./62.8), apert=CircularAperture(lc1, 12.7))
cementsurf = Surface(lc2, shape=surfShape.Conic(lc2, curv=-1./45.7), apert=CircularAperture(lc2, 12.7))
rearsurf = Surface(lc3, shape=surfShape.Conic(lc3, curv=-1./128.2), apert=CircularAperture(lc3, 12.7))
image = Surface(lc4)


elem = OpticalElement(lc0, name="thorlabs_AC_254-100-A")

bk7 = ConstantIndexGlass(lc1, n=1.5168)
sf5 = ConstantIndexGlass(lc2, n=1.6727)

elem.addMaterial("BK7", bk7)
elem.addMaterial("SF5", sf5)

elem.addSurface("stop", stopsurf, (None, None))
elem.addSurface("front", frontsurf, (None, "BK7"))
elem.addSurface("cement", cementsurf, ("BK7", "SF5"))
elem.addSurface("rear", rearsurf, ("SF5", None))
elem.addSurface("image", image, (None, None))

s.addElement("AC254-100", elem)

rstobj = raster.MeridionalFan()
(px, py) = rstobj.getGrid(20)

rpup = 11.43
o = np.vstack((rpup*px, rpup*py, -5.*np.ones_like(px)))

k = np.zeros_like(o)
k[2,:] = 1.0 #2.*math.pi/wavelength

ey = np.zeros_like(o)
ey[1,:] =  1.

E0 = np.cross(k, ey, axisa=0, axisb=0).T

sysseq = [("AC254-100", [("stop", {"is_stop":True}), ("front", {}), ("cement", {}), ("rear", {}), ("image", {})])]

phi = 5.*math.pi/180.0

initialbundle = RayBundle(x0=o, k0=k, Efield0=E0, wave=wavelength)
r2 = s.seqtrace(initialbundle, sysseq)

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

for r in r2:
    r.draw2d(ax, color="blue", plane_normal=pn, up=up) 

s.draw2d(ax, color="grey", vertices=50, plane_normal=pn, up=up) # try for phi=0.
#s.draw2d(ax, color="grey", inyzplane=False, vertices=50, plane_normal=pn, up=up) # try for phi=pi/4


plt.show()


