#!/usr/bin/env/python
"""
Pyrate - Optical raytracing based on Python

Copyright (C) 2014-2020
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
import math
import logging

from pyrateoptics.sampling2d import raster
from pyrateoptics.raytracer.material.material_grin import IsotropicGrinMaterial
from pyrateoptics.raytracer.surface_shape import Conic
from pyrateoptics.raytracer.optical_element import OpticalElement
from pyrateoptics.raytracer.surface import Surface
from pyrateoptics.raytracer.optical_system import OpticalSystem

from pyrateoptics.raytracer.aperture import CircularAperture
from pyrateoptics.raytracer.localcoordinates import LocalCoordinates

from pyrateoptics import raytrace, draw

from pyrateoptics.core.serializer import Serializer, Deserializer

from pprint import pprint

logging.basicConfig(level=logging.DEBUG)

wavelength = 0.5876e-3

# definition of optical system
s = OpticalSystem.p()

lc0 = s.addLocalCoordinateSystem(
    LocalCoordinates.p(name="obj", decz=0.0),
    refname=s.rootcoordinatesystem.name)
lc1 = s.addLocalCoordinateSystem(
    LocalCoordinates.p(name="surf1", decz=10.0, tiltx=5.*math.pi/180.0),
    refname=lc0.name)  # objectDist
lc2 = s.addLocalCoordinateSystem(
    LocalCoordinates.p(name="surf2", decz=20.0, tiltx=10.*math.pi/180.0),
    refname=lc1.name)
lc3 = s.addLocalCoordinateSystem(
    LocalCoordinates.p(name="image", decz=10.0), refname=lc2.name)


stopsurf = Surface.p(lc0)
surf1 = Surface.p(lc1, shape=Conic.p(lc1, curv=1./24.0),
                  aperture=CircularAperture.p(lc1, maxradius=5.0))
surf2 = Surface.p(lc2, shape=Conic.p(lc2, curv=-1./24.0),
                  aperture=CircularAperture.p(lc2, maxradius=5.0))
image = Surface.p(lc3)

elem = OpticalElement.p(lc0, name="grinelement")

mysource =\
r"""

import numpy as np

grin_strength = 0.5


def nfunc(x, **kw):
    '''
    Refractive index function.
    '''
    return grin_strength*np.exp(-x[0]**2 - 4.*x[1]**2)+1.0
    # (2.5 - (x**2 + 100.0*y**4)/10.**2)


def dndx(x, **kw):
    '''
    d/dx of refractive index function
    '''
    return -2.*x[0]*grin_strength*np.exp(-x[0]**2 - 4.*x[1]**2)
    # -2*x/10.**2


def dndy(x, **kw):
    '''
    d/dy of refractive index function
    '''
    return -2.*4.*x[1]*grin_strength*np.exp(-x[0]**2 - 4.*x[1]**2)
    # -100.0*4.0*y**3/10.**2


def dndz(x, **kw):
    '''
    d/dz of refractive index function
    '''
    return np.zeros_like(x[0])


def bnd(x):
    '''
    Boundary function.
    '''
    return x[0]**2 + x[1]**2 < 10.**2
"""

# grinmaterial = ConstantIndexGlass.p(lc1, 1.0 + grin_strength)
grinmaterial = IsotropicGrinMaterial.p(lc1, mysource,
                                       "nfunc",
                                       "dndx",
                                       "dndy",
                                       "dndz",
                                       "bnd",
                                       parameterlist=[("n0", 0.5)])
grinmaterial.annotations["ds"] = 0.05
grinmaterial.annotations["energyviolation"] = 0.01

elem.addMaterial("grin", grinmaterial)

elem.addSurface("object", stopsurf, (None, None))
elem.addSurface("surf1", surf1, (None, "grin"))
elem.addSurface("surf2", surf2, ("grin", None))
elem.addSurface("image", image, (None, None))

s.addElement("grinelement", elem)

sysseq = [("grinelement",
           [("object", {"is_stop": True}),
            ("surf1", {}), ("surf2", {}),
            ("image", {})])]


r2 = raytrace(s, sysseq, 21,
              {"startz": -5., "radius": 2.5, "raster": raster.MeridionalFan()},
              wave=wavelength)
draw(s, r2)
