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

from core import raster
from core import material
from core import surfShape
from core.optical_element import OpticalSystemNew, SurfaceNew, OpticalElement
from core.ray import RayBundleNew

from core.aperture import CircularAperture
from core.coordinates import LocalCoordinates

from core.globalconstants import canonical_ex, canonical_ey

import math

wavelength = 0.5876e-3

# definition of optical system
s = OpticalSystemNew() 

lc0 = s.addLocalCoordinateSystem(LocalCoordinates(name="stop", decz=0.0), refname=s.rootcoordinatesystem.name)
lc1 = s.addLocalCoordinateSystem(LocalCoordinates(name="surf1", decz=-1.048), refname=lc0.name) # objectDist
lc2 = s.addLocalCoordinateSystem(LocalCoordinates(name="surf2", decz=4.0), refname=lc1.name)
lc3 = s.addLocalCoordinateSystem(LocalCoordinates(name="surf3", decz=2.5), refname=lc2.name)
lc4 = s.addLocalCoordinateSystem(LocalCoordinates(name="image", decz=97.2), refname=lc3.name)


stopsurf = SurfaceNew(lc0)
frontsurf = SurfaceNew(lc1, shape=surfShape.Conic(lc1, curv=1./62.8), apert=CircularAperture(lc1, 12.7))
cementsurf = SurfaceNew(lc2, shape=surfShape.Conic(lc2, curv=-1./45.7), apert=CircularAperture(lc2, 12.7))
rearsurf = SurfaceNew(lc3, shape=surfShape.Conic(lc3, curv=-1./128.2), apert=CircularAperture(lc3, 12.7))
image = SurfaceNew(lc4)


elem = OpticalElement(lc0, label="thorlabs_AC_254-100-A")

bk7 = material.ConstantIndexGlass(lc1, n=1.5168)
sf5 = material.ConstantIndexGlass(lc2, n=1.6727)

elem.addMaterial("BK7", bk7)
elem.addMaterial("SF5", sf5)

elem.addSurface("stop", stopsurf, (None, None))
elem.addSurface("front", frontsurf, (None, "BK7"))
elem.addSurface("cement", cementsurf, ("BK7", "SF5"))
elem.addSurface("rear", rearsurf, ("SF5", None))
elem.addSurface("image", image, (None, None))

s.addElement("AC254-100", elem)

print(s.rootcoordinatesystem.pprint())

rstobj = raster.RectGrid()
(px, py) = rstobj.getGrid(100)

rpup = 11.43
o = np.vstack((rpup*px, rpup*py, -5.*np.ones_like(px)))

k = np.zeros_like(o)
k[2,:] = 2.*math.pi/wavelength

ey = np.zeros_like(o)
ey[1,:] =  1.

E0 = np.cross(k, ey, axisa=0, axisb=0).T

sysseq = [("AC254-100", ["stop", "front", "cement", "rear", "image"])]

phi = 5.*math.pi/180.0

initialbundle = RayBundleNew(x0=o, k0=k, Efield0=E0, wave=wavelength)
r2 = s.seqtrace(initialbundle, sysseq)

pilotbundle = RayBundleNew(
                x0 = np.array([[0], [0], [0]]), 
                k0 = np.array([[0], [2.*math.pi/wavelength*math.sin(phi)], [2.*math.pi/wavelength*math.cos(phi)]]), 
                Efield0 = np.array([[1], [0], [0]]), wave=wavelength
                )

pilotray = s.seqtrace(pilotbundle, sysseq)

#for (ind, r) in enumerate(pilotray.raybundles):
#    print("pilot bundle %d" % (ind,))
#    print(r.x)
    
#print("last coordinates")
#print(r2.raybundles[-1].x[-1, :, :])


fig = plt.figure(1)
ax = fig.add_subplot(111)

ax.axis('equal')
ax.set_axis_bgcolor('black')

phi = math.pi/4
pn = np.array([math.cos(phi), 0, math.sin(phi)]) # canonical_ex
up = canonical_ey

print("drawing!")
r2.draw2d(ax, color="blue", plane_normal=pn, up=up) 
pilotray.draw2d(ax, color="green", plane_normal=pn, up=up)
for e in s.elements.itervalues():
    for surfs in e.surfaces.itervalues():
        #surfs.draw2d(ax, color="grey", vertices=50, plane_normal=pn, up=up) # try for phi=0.
        surfs.draw2d(ax, color="grey", inyzplane=False, vertices=50, plane_normal=pn, up=up) # try for phi=pi/4

#plots.drawLayout2d(ax, s, [pilotpath])
#plots.drawLayout2d(ax, s, [r2])

plt.show()


