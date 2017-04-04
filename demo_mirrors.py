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
lc1 = s.addLocalCoordinateSystem(LocalCoordinates(name="m1", decz=50.0, tiltx=-math.pi/8), refname=lc0.name) # objectDist
lc2 = s.addLocalCoordinateSystem(LocalCoordinates(name="m2", decz=-50.0, tiltx=math.pi/8), refname=lc1.name)
lc3 = s.addLocalCoordinateSystem(LocalCoordinates(name="m3", decz=50.0, decy=-20), refname=lc2.name)
lc4 = s.addLocalCoordinateSystem(LocalCoordinates(name="image", decz=-50, decy=-20), refname=lc3.name)


stopsurf = SurfaceNew(lc0)
frontsurf = SurfaceNew(lc1, shape=surfShape.Conic(lc1), apert=CircularAperture(lc1, 12.7))
cementsurf = SurfaceNew(lc2, shape=surfShape.Conic(lc2), apert=CircularAperture(lc2, 12.7))
rearsurf = SurfaceNew(lc3, shape=surfShape.Conic(lc3), apert=CircularAperture(lc3, 12.7))
image = SurfaceNew(lc4)


elem = OpticalElement(lc0, label="TMA")

elem.addSurface("stop", stopsurf, (None, None))
elem.addSurface("m1", frontsurf, (None, None))
elem.addSurface("m2", cementsurf, (None, None))
elem.addSurface("m3", rearsurf, (None, None))
elem.addSurface("image", image, (None, None))

s.addElement("TMA", elem)

print(s.rootcoordinatesystem.pprint())

rstobj = raster.MeridionalFan()
(px, py) = rstobj.getGrid(5)

rpup = 10
o = np.vstack((rpup*px, rpup*py, -5.*np.ones_like(px)))

k = np.zeros_like(o)
k[2,:] = 2.*math.pi/wavelength

ey = np.zeros_like(o)
ey[1,:] =  1.

E0 = np.cross(k, ey, axisa=0, axisb=0).T

sysseq = [("TMA", [("stop", True, True), ("m1", False, True), ("m2", False, True), ("m3", False, True), ("image", True, True)])] 

phi = 5.*math.pi/180.0

initialbundle = RayBundleNew(x0=o, k0=k, Efield0=E0, wave=wavelength)
r2 = s.seqtrace(initialbundle, sysseq)

#pilotbundle = RayBundleNew(
#                x0 = np.array([[0], [0], [0]]), 
#                k0 = np.array([[0], [2.*math.pi/wavelength*math.sin(phi)], [2.*math.pi/wavelength*math.cos(phi)]]), 
#                Efield0 = np.array([[1], [0], [0]]), wave=wavelength
#                )

#pilotray = s.seqtrace(pilotbundle, sysseq)

for (ind, r) in enumerate(r2.raybundles):
    print("bundle %d" % (ind,))
    print(r.x)
    print(r.k)
    
#print("last coordinates")
#print(r2.raybundles[-1].x[-1, :, :])


fig = plt.figure(1)
ax = fig.add_subplot(111)

ax.axis('equal')
ax.set_axis_bgcolor('black')

phi = 0. #math.pi/4
pn = np.array([math.cos(phi), 0, math.sin(phi)]) # canonical_ex
up = canonical_ey

print("drawing!")
r2.draw2d(ax, color="blue", plane_normal=pn, up=up) 
#pilotray.draw2d(ax, color="green", plane_normal=pn, up=up)
for e in s.elements.itervalues():
    for surfs in e.surfaces.itervalues():
        surfs.draw2d(ax, color="grey", vertices=50, plane_normal=pn, up=up) # try for phi=0.
        #surfs.draw2d(ax, color="grey", inyzplane=False, vertices=50, plane_normal=pn, up=up) # try for phi=pi/4

#plots.drawLayout2d(ax, s, [pilotpath])
#plots.drawLayout2d(ax, s, [r2])

plt.show()

