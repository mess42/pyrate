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
from pyrateoptics.material.material_isotropic import ModelGlass
from pyrateoptics.raytracer import surfShape
from pyrateoptics.raytracer.optical_element import OpticalElement
from pyrateoptics.raytracer.surface import Surface
from pyrateoptics.raytracer.optical_system import OpticalSystem
from pyrateoptics.raytracer.ray import RayBundle

from pyrateoptics.raytracer.aperture import CircularAperture
from pyrateoptics.raytracer.localcoordinates import LocalCoordinates

from pyrateoptics import collimated_bundle, draw

import math
import logging
logging.basicConfig(level=logging.DEBUG)

wavelength = 0.5876e-3

wave_red = 0.700e-3
wave_blue = 0.470e-3

# definition of optical system
s = OpticalSystem() 

deg = math.pi/180.

lc0 = s.addLocalCoordinateSystem(LocalCoordinates(name="stop", decz=0.0), refname=s.rootcoordinatesystem.name)
lccomprism = s.addLocalCoordinateSystem(LocalCoordinates(name="prismcenter", decz=50.0), refname=lc0.name)

lc1 = s.addLocalCoordinateSystem(LocalCoordinates(name="surf1", decz=-10.0, tiltx=30.*deg), refname=lccomprism.name) # objectDist
lc2 = s.addLocalCoordinateSystem(LocalCoordinates(name="surf2", decz=10.0, tiltx=-30.*deg), refname=lccomprism.name)
lc3 = s.addLocalCoordinateSystem(LocalCoordinates(name="image", decz=50.0), refname=lccomprism.name)


stopsurf = Surface(lc0)
frontsurf = Surface(lc1, shape=surfShape.Conic(lc1, curv=0), apert=CircularAperture(lc1, 20.0))
rearsurf = Surface(lc2, shape=surfShape.Conic(lc2, curv=0), apert=CircularAperture(lc2, 20.0))
image = Surface(lc3)


elem = OpticalElement(lc0, name="prism")

glass = ModelGlass(lc1)


elem.addMaterial("glass", glass)

elem.addSurface("stop", stopsurf, (None, None))
elem.addSurface("surf1", frontsurf, (None, "glass"))
elem.addSurface("surf2", rearsurf, ("glass", None))
elem.addSurface("image", image, (None, None))

s.addElement("prism", elem)

sysseq = [("prism", 
               [("stop", {"is_stop":True}), 
                ("surf1", {}), 
                ("surf2", {}), 
                ("image", {})])]

(o, k_red, E0_red) = collimated_bundle(20, {"opticalsystem": s, "radius": 5.0, "startz":-5., "starty":-20., "anglex":23*deg, "raster":raster.MeridionalFan()}, wave=wave_red)
(o, k_blue, E0_blue) = collimated_bundle(20, {"opticalsystem": s, "radius": 5.0, "startz":-5., "starty":-20., "anglex":23*deg, "raster":raster.MeridionalFan()}, wave=wave_blue)

initialbundle_red = RayBundle(x0=o, k0=k_red, Efield0=E0_red, wave=wave_red)
initialbundle_blue = RayBundle(x0=o, k0=k_blue, Efield0=E0_blue, wave=wave_blue)
r_red = s.seqtrace(initialbundle_red, sysseq)
r_blue = s.seqtrace(initialbundle_blue, sysseq)

draw(s, [(r_red, "red"), (r_blue, "blue")])

