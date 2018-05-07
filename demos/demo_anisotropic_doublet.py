#!/usr/bin/env/python
"""
Pyrate - Optical raytracing based on Python

Copyright (C) 2014-2018
               by     Moritz Esslinger moritz.esslinger@web.de
               and    Johannes Hartung j.hartung@gmx.net
               and    Uwe Lippmann  uwe.lippmann@web.de
               and    Thomas Heinze t.heinze@uni-jena.de
               and    others

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


from pyrateoptics.sampling2d import raster
from pyrateoptics.material.material_anisotropic import AnisotropicMaterial
from pyrateoptics.raytracer import surfShape
from pyrateoptics.raytracer.optical_element import OpticalElement
from pyrateoptics.raytracer.surface import Surface
from pyrateoptics.raytracer.optical_system import OpticalSystem
from pyrateoptics.raytracer.ray import RayBundle

from pyrateoptics.raytracer.aperture import CircularAperture
from pyrateoptics.raytracer.localcoordinates import LocalCoordinates

from pyrateoptics import draw, collimated_bundle

import logging
logging.basicConfig(level=logging.DEBUG)

wavelength = 0.5876e-3

# definition of optical system
s = OpticalSystem(name='os') 

lc0 = s.addLocalCoordinateSystem(LocalCoordinates(name="stop", decz=0.0), refname=s.rootcoordinatesystem.name)
lc1 = s.addLocalCoordinateSystem(LocalCoordinates(name="surf1", decz=-1.048), refname=lc0.name) # objectDist
lc2 = s.addLocalCoordinateSystem(LocalCoordinates(name="surf2", decz=4.0), refname=lc1.name)
lc3 = s.addLocalCoordinateSystem(LocalCoordinates(name="surf3", decz=2.5), refname=lc2.name)
lc4 = s.addLocalCoordinateSystem(LocalCoordinates(name="image", decz=97.2), refname=lc3.name)


stopsurf = Surface(lc0, name="stopsurf")
frontsurf = Surface(lc1, name="frontsurf", shape=surfShape.Conic(lc1, curv=1./62.8, name='conic1'), apert=CircularAperture(lc1, 12.7))
cementsurf = Surface(lc2, name="cementsurf", shape=surfShape.Conic(lc2, curv=-1./45.7, name='conic2'), apert=CircularAperture(lc2, 12.7))
rearsurf = Surface(lc3, name="rearsurf", shape=surfShape.Conic(lc3, curv=-1./128.2, name='conic3'), apert=CircularAperture(lc3, 12.7))
image = Surface(lc4, name="imagesurf")


elem = OpticalElement(lc0, name="thorlabs_AC_254-100-A")

rnd_data1 = np.random.random((3, 3)) #np.eye(3)
rnd_data2 = np.random.random((3, 3))#np.zeros((3, 3))#
rnd_data3 = np.random.random((3, 3)) #np.eye(3)
rnd_data4 = np.random.random((3, 3))#np.zeros((3, 3))#

# isotropic tests

#bk7 = material.ConstantIndexGlass(lc1, n=1.5168)
#sf5 = material.ConstantIndexGlass(lc2, n=1.6727)

myeps1 = 1.5168**2*np.eye(3)
myeps2 = 1.6727**2*np.eye(3)

# anisotropic materials

#myeps1 = rnd_data1 + complex(0, 1)*rnd_data2
#myeps2 = rnd_data3 + complex(0, 1)*rnd_data4

crystal1 = AnisotropicMaterial(lc1, myeps1, name="crystal1")
crystal2 = AnisotropicMaterial(lc2, myeps2, name="crystal2")


elem.addMaterial("crystal1", crystal1)
elem.addMaterial("crystal2", crystal2)

elem.addSurface("stop", stopsurf, (None, None))
elem.addSurface("front", frontsurf, (None, "crystal1"))
elem.addSurface("cement", cementsurf, ("crystal1", "crystal2"))
elem.addSurface("rear", rearsurf, ("crystal2", None))
elem.addSurface("image", image, (None, None))

s.addElement("AC254-100", elem)

sysseq = [("AC254-100", [("stop", {}), ("front", {}), ("cement", {}), ("rear", {}), ("image", {})])]


(o, k, E0) = collimated_bundle(10, {"opticalsystem":s, "radius":11.43, "startz":-5., "raster":raster.MeridionalFan()}, wave=wavelength)
initialbundle = RayBundle(x0=o, k0=k, Efield0=E0, wave=wavelength)
r2 = s.seqtrace(initialbundle, sysseq, splitup=True)

draw(s, [(r2[0], "blue"), (r2[1], "green")])

