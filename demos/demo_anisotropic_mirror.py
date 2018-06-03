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
import math
import logging
logging.basicConfig(level=logging.DEBUG)


from pyrateoptics.sampling2d import raster
from pyrateoptics.material.material_anisotropic import AnisotropicMaterial
from pyrateoptics.raytracer import surfShape
from pyrateoptics.raytracer.optical_element import OpticalElement
from pyrateoptics.raytracer.surface import Surface
from pyrateoptics.raytracer.optical_system import OpticalSystem
from pyrateoptics.raytracer.ray import RayBundle

from pyrateoptics.raytracer.aperture import CircularAperture
from pyrateoptics.raytracer.localcoordinates import LocalCoordinates
from pyrateoptics.raytracer.globalconstants import degree


from pyrateoptics import draw, raytrace



wavelength = 0.5876e-3

# definition of optical system
s = OpticalSystem(name='os') 

lc0 = s.addLocalCoordinateSystem(LocalCoordinates(name="stop", decz=1.0), refname=s.rootcoordinatesystem.name)
lc1 = s.addLocalCoordinateSystem(LocalCoordinates(name="surf1", decz=10.0), refname=lc0.name) # objectDist
lc2 = s.addLocalCoordinateSystem(LocalCoordinates(name="surf2", decz=5.0, tiltx=10*math.pi/180.0), refname=lc1.name)
lc3 = s.addLocalCoordinateSystem(LocalCoordinates(name="image", decz=-5.0, tiltx=-10*math.pi/180.0), refname=lc2.name)


stopsurf = Surface(lc0)
frontsurf = Surface(lc1, shape=surfShape.Conic(lc1, curv=0), apert=CircularAperture(lc1, 10.0))
rearsurf = Surface(lc2, shape=surfShape.Conic(lc2, curv=0), apert=CircularAperture(lc3, 10.0))
image = Surface(lc3)


elem = OpticalElement(lc0, name="crystalelem")

no = 1.5
neo = 1.8

myeps = np.array([[no, 0, 0], [0, no, 0], [0, 0, neo]]) 

crystal = AnisotropicMaterial(lc1, myeps)


elem.addMaterial("crystal", crystal)

elem.addSurface("stop", stopsurf, (None, None))
elem.addSurface("front", frontsurf, (None, "crystal"))
elem.addSurface("rear", rearsurf, ("crystal", "crystal"))
elem.addSurface("image", image, ("crystal", None))

s.addElement("crystalelem", elem)

sysseq = [
    ("crystalelem", [
        ("stop", {}), 
        ("front", {}), 
        ("rear", {"is_mirror":True}), 
        ("image", {})]
    )]

rays = raytrace(s, sysseq, 10,\
         {"radius":20*degree, "startz":-5., "raster":raster.MeridionalFan()},\
         bundletype="divergent", traceoptions={"splitup":True}, wave=wavelength)


draw(s, rays)

