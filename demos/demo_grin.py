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

import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from distutils.version import StrictVersion

from core import raster
from core.material_grin import IsotropicGrinMaterial
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

lc0 = s.addLocalCoordinateSystem(LocalCoordinates(name="obj", decz=0.0), refname=s.rootcoordinatesystem.name)
lc1 = s.addLocalCoordinateSystem(LocalCoordinates(name="surf1", decz=10.0, tiltx=5.*math.pi/180.0), refname=lc0.name) # objectDist
lc2 = s.addLocalCoordinateSystem(LocalCoordinates(name="surf2", decz=20.0, tiltx=10.*math.pi/180.0), refname=lc1.name)
lc3 = s.addLocalCoordinateSystem(LocalCoordinates(name="image", decz=10.0), refname=lc2.name)


stopsurf = Surface(lc0)
surf1 = Surface(lc1, shape=surfShape.Conic(lc1, curv=1./24.0), apert=CircularAperture(lc1, 5.0))
surf2 = Surface(lc2, shape=surfShape.Conic(lc2, curv=-1./24.0), apert=CircularAperture(lc2, 5.0))
image = Surface(lc3)

elem = OpticalElement(lc0, name="grinelement")

grin_strength = 0.5


def nfunc(x):
    return grin_strength*np.exp(-x[0]**2 - 4.*x[1]**2)+1.0#(2.5 - (x**2 + 100.0*y**4)/10.**2)

def dndx(x):
    return -2.*x[0]*grin_strength*np.exp(-x[0]**2 - 4.*x[1]**2)#-2*x/10.**2

def dndy(x):
    return -2.*4.*x[1]*grin_strength*np.exp(-x[0]**2 - 4.*x[1]**2) #-100.0*4.0*y**3/10.**2

def dndz(x):
    return np.zeros_like(x[0])

def bnd(x):
    return x[0]**2 + x[1]**2 < 10.**2

#grinmaterial = ConstantIndexGlass(lc1, 1.0 + grin_strength) 
grinmaterial = IsotropicGrinMaterial(lc1, nfunc, dndx, dndy, dndz, bnd, ds=0.05, energyviolation=0.01)

elem.addMaterial("grin", grinmaterial)

elem.addSurface("object", stopsurf, (None, None))
elem.addSurface("surf1", surf1, (None, "grin"))
elem.addSurface("surf2", surf2, ("grin", None))
elem.addSurface("image", image, (None, None))

s.addElement("grinelement", elem)

rstobj = raster.MeridionalFan()
(px, py) = rstobj.getGrid(21)

rpup = 2.5
o = np.vstack((rpup*px, rpup*py, -5.*np.ones_like(px)))

k = np.zeros_like(o)
k[2,:] = 1. #2.*math.pi/wavelength

ey = np.zeros_like(o)
ey[1,:] =  1.

E0 = np.cross(k, ey, axisa=0, axisb=0).T

sysseq = [("grinelement", [("object", {"is_stop":True}), ("surf1", {}), ("surf2", {}), ("image", {})])]

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

for r in r2:
    r.draw2d(ax, color="blue", plane_normal=pn, up=up)
s.draw2d(ax, color="grey", vertices=50, plane_normal=pn, up=up) # try for phi=0.
#s.draw2d(ax, color="grey", inyzplane=False, vertices=50, plane_normal=pn, up=up) # try for phi=pi/4


plt.show()

