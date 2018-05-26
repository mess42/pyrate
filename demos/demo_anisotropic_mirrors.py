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
from pyrateoptics.raytracer.optical_system import OpticalSystem
from pyrateoptics.raytracer.surface import Surface
from pyrateoptics.raytracer.ray import RayBundle

from pyrateoptics.raytracer.aperture import CircularAperture
from pyrateoptics.raytracer.localcoordinates import LocalCoordinates

from pyrateoptics import draw, collimated_bundle


import pyrateoptics.raytracer.helpers

wavelength = 0.5876e-3

rnd_data1 = np.random.random((3, 3)) #np.eye(3)
rnd_data2 = np.random.random((3, 3))#np.zeros((3, 3))#
lc = LocalCoordinates("1")
myeps = np.eye(3) + 0.1*rnd_data1 + 0.01*complex(0, 1)*rnd_data2 # aggressive complex choice of myeps
#myeps = np.eye(3) + 0.01*np.random.random((3, 3))
crystal = AnisotropicMaterial(lc, myeps)


# definition of optical system
s = OpticalSystem(matbackground=crystal) 

lc0 = s.addLocalCoordinateSystem(LocalCoordinates(name="object", decz=0.0), refname=s.rootcoordinatesystem.name)
lc1 = s.addLocalCoordinateSystem(LocalCoordinates(name="m1", decz=50.0, tiltx=-math.pi/8), refname=lc0.name) # objectDist
lc2 = s.addLocalCoordinateSystem(LocalCoordinates(name="m2_stop", decz=-50.0, decy=-20, tiltx=math.pi/16), refname=lc1.name)
lc3 = s.addLocalCoordinateSystem(LocalCoordinates(name="m3", decz=50.0, decy=-30, tiltx=3*math.pi/32), refname=lc2.name)
lc4 = s.addLocalCoordinateSystem(LocalCoordinates(name="image1", decz=-50, decy=-15, tiltx=-math.pi/16), refname=lc3.name)
lc5 = s.addLocalCoordinateSystem(LocalCoordinates(name="oapara", decz=-100, decy=-35), refname=lc4.name)
lc5ap = s.addLocalCoordinateSystem(LocalCoordinates(name="oaparaap", decz=0, decy=35), refname=lc5.name)
lc6 = s.addLocalCoordinateSystem(LocalCoordinates(name="image2", decz=55, tiltx=1*math.pi/32), refname=lc5.name)

objectsurf = Surface(lc0)
m1surf = Surface(lc1, shape=surfShape.Conic(lc1, curv=-0.01), apert=CircularAperture(lc1, 20.))
m2surf = Surface(lc2, shape=surfShape.Conic(lc2, curv=0.01), apert=CircularAperture(lc2, 12.7))
m3surf = Surface(lc3, shape=surfShape.Conic(lc3, curv=-0.006), apert=CircularAperture(lc3, 20.7))
image1 = Surface(lc4)
oapara = Surface(lc3, shape=surfShape.Conic(lc5, curv=0.01, cc=-1.), apert=CircularAperture(lc5ap, 30.0))
image2 = Surface(lc6, apert=CircularAperture(lc6, 20.0))


elem = OpticalElement(lc0, name="TMA")

#elem.addMaterial("crystal", crystal)

elem.addSurface("object", objectsurf, (None, None))
elem.addSurface("m1", m1surf, (None, None))
elem.addSurface("m2", m2surf, (None, None))
elem.addSurface("m3", m3surf, (None, None))
elem.addSurface("image1", image1, (None, None))
elem.addSurface("oapara", oapara, (None, None))
elem.addSurface("image2", image2, (None, None))

s.addElement("TMA", elem)

print(s.rootcoordinatesystem.pprint())

(o, k, E0) = collimated_bundle(5, {"opticalsystem":s, "startz":-5., "radius":10., "raster":raster.MeridionalFan()}, wave=wavelength)

sysseq = [("TMA", [
            ("object", {}), 
            ("m1", {"is_mirror":True}), 
            ("m2", {"is_mirror":True}), 
            ("m3", {"is_mirror":True}), 
            ("image1", {}), 
            ("oapara", {"is_mirror":True}), 
            ("image2", {}) ])] 

sysseq_pilot = [("TMA", 
                 [
                    ("object", {}), 
                    ("m1", {"is_mirror":True}), 
                    ("m2", {"is_mirror":True}), 
                    ("m3", {"is_mirror":True}), 
                    ("m2", {"is_mirror":True}),
                    ("m1", {"is_mirror":True}),
                    ("m2", {"is_mirror":True}),
                    ("m1", {"is_mirror":True}),
                    ("m2", {"is_mirror":True})
                ])
                ] 
                
obj_dx = 0.1
obj_dphi = 1.*math.pi/180.0

kwave = 2.*math.pi/wavelength

initialbundle = RayBundle(x0=o, k0=k, Efield0=E0, wave=wavelength)
r2 = s.seqtrace(initialbundle, sysseq)

#pilotbundle = RayBundle(
#                x0 = np.array([[0], [0], [0]]), 
#                k0 = np.array([[0], [kwave*math.sin(phi)], [kwave*math.cos(phi)]]), 
#                Efield0 = np.array([[1], [0], [0]]), wave=wavelength
#                )
#pilotray = s.seqtrace(pilotbundle, sysseq_pilot)

kw = 5*math.pi/180.

pilotbundles = pyrateoptics.raytracer.helpers.build_pilotbundle(objectsurf, crystal, (obj_dx, obj_dx), (obj_dphi, obj_dphi), kunitvector=np.array([0, math.sin(kw), math.cos(kw)]), num_sampling_points=3)
(pilotray2, r3) = s.para_seqtrace(pilotbundles[-1], initialbundle, sysseq)

draw(s, [(r2, "blue"), (r3, "orange"), (pilotray2, "red")])
