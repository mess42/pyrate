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
from core.optical_element import OpticalSystemNew, SurfaceNew, OpticalElement
from core.ray import RayPath, RayBundle, RayBundleNew

from core import plots
from core.aperture import CircularAperture
from core.coordinates import LocalCoordinates

import math

wavelength = 0.5876e-3

# definition of optical system
s = OpticalSystemNew() 

lc0 = s.addLocalCoordinateSystem(LocalCoordinates(name="stop", decz=0.0))
lc1 = s.addLocalCoordinateSystem(LocalCoordinates(name="surf1", decz=-1.048)) # objectDist
lc2 = s.addLocalCoordinateSystem(LocalCoordinates(name="surf2", decz=4.0))
lc3 = s.addLocalCoordinateSystem(LocalCoordinates(name="surf3", decz=2.5))
lc4 = s.addLocalCoordinateSystem(LocalCoordinates(name="image", decz=97.2))


stopsurf = SurfaceNew(lc0)
frontsurf = SurfaceNew(lc1, surfShape.Conic(lc1, curv=1./62.8), aperture=CircularAperture(lc1, 12.7))
cementsurf = SurfaceNew(lc2, surfShape.Conic(lc2, curv=-1./45.7), aperture=CircularAperture(lc2, 12.7))
rearsurf = SurfaceNew(lc3, surfShape.Conic(lc3, curv=-1./128.2), aperture=CircularAperture(lc3, 12.7))
image = SurfaceNew(lc4)


elem = OpticalElement(lc0, label="thorlabs_AC_254-100-A")

elem.addMaterial("BK7", material.ConstantIndexGlass(lc1, n=1.5168))
elem.addMaterial("SF5", material.ConstantIndexGlass(lc2, n=1.6727))

elem.addSurface("stop", stopsurf, (None, None))
elem.addSurface("front", frontsurf, (None, "BK7"))
elem.addSurface("cement", cementsurf, ("BK7", "SF5"))
elem.addSurface("rear", rearsurf, ("SF5", None))
elem.addSurface("image", image, (None, None))

s.addElement("AC254-100", elem)


rstobj = raster.RectGrid()
(px, py) = rstobj.getGrid(10)

rpup = 11.43
o = np.vstack((rpup*px, rpup*py, -5.*np.ones_like(px)))

k = np.zeros_like(o)
k[2,:] = 2.*math.pi/wavelength

ey = np.zeros_like(o)
ey[1,:] =  1.

E0 = np.cross(k, ey, axisa=0, axisb=0).T


initialbundle = RayBundleNew(x0=o, k0=k, Efield0=E0, wave=wavelength)

r2 = s.seqtrace(initialbundle, [("AC254-100", ["stop", "front", "cement", "rear", "image"])])

for (ind, r) in enumerate(r2.raybundles):
    print("bundle %d" % (ind,))
    print(r.x)
    print(r.valid)
    #print(r.k)
    
print(r2.raybundles[-1].x[-1, :, :])

"""
fig = plt.figure(1)
ax = fig.add_subplot(111)

ax.axis('equal')
ax.set_axis_bgcolor('black')

#plots.drawLayout2d(ax, s, [pilotpath])
plots.drawLayout2d(ax, s, [r2])

print()

plt.show()
"""

