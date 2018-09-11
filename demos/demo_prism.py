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

import logging


from pyrateoptics.sampling2d import raster
from pyrateoptics.material.material_isotropic import ModelGlass
from pyrateoptics.raytracer.surface_shape import Conic
from pyrateoptics.raytracer.optical_element import OpticalElement
from pyrateoptics.raytracer.surface import Surface
from pyrateoptics.raytracer.optical_system import OpticalSystem

from pyrateoptics.raytracer.aperture import CircularAperture
from pyrateoptics.raytracer.localcoordinates import LocalCoordinates

from pyrateoptics import raytrace, draw
from pyrateoptics.raytracer.globalconstants import degree

logging.basicConfig(level=logging.DEBUG)

wavelength = 0.5876e-3

wave_red = 0.700e-3
wave_blue = 0.470e-3

# definition of optical system
s = OpticalSystem()

lc0 = s.addLocalCoordinateSystem(
            LocalCoordinates(name="stop", decz=0.0),
            refname=s.rootcoordinatesystem.name)
lccomprism = s.addLocalCoordinateSystem(
            LocalCoordinates(name="prismcenter", decz=50.0),
            refname=lc0.name)

lc1 = s.addLocalCoordinateSystem(
            LocalCoordinates(name="surf1", decz=-10.0, tiltx=30.*degree),
            refname=lccomprism.name)  # objectDist
lc2 = s.addLocalCoordinateSystem(
            LocalCoordinates(name="surf2", decz=10.0, tiltx=-30.*degree),
            refname=lccomprism.name)
lc3 = s.addLocalCoordinateSystem(
            LocalCoordinates(name="image", decz=50.0),
            refname=lccomprism.name)


stopsurf = Surface(lc0)
frontsurf = Surface(lc1, shape=Conic(lc1, curv=0),
                    aperture=CircularAperture(lc1, maxradius=20.0))
rearsurf = Surface(lc2, shape=Conic(lc2, curv=0),
                   aperture=CircularAperture(lc2, maxradius=20.0))
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
          [("stop", {"is_stop": True}),
           ("surf1", {}),
           ("surf2", {}),
           ("image", {})])]

raysdict = {"radius": 5.0, "startz": -5., "starty": -20., "anglex": 23*degree,
            "raster": raster.MeridionalFan()}

r_red = raytrace(s, sysseq, 20, raysdict, wave=wave_red)[0]
r_blue = raytrace(s, sysseq, 20, raysdict, wave=wave_blue)[0]


draw(s, [(r_red, "red"), (r_blue, "blue")])
