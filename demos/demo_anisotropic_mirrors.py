#!/usr/bin/env python3
# -*- coding: utf-8 -*-
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

import numpy as np

from pyrateoptics.sampling2d import raster
from pyrateoptics.raytracer.material.material_anisotropic import\
     AnisotropicMaterial
from pyrateoptics.raytracer.surface_shape import Conic
from pyrateoptics.raytracer.optical_element import OpticalElement
from pyrateoptics.raytracer.optical_system import OpticalSystem
from pyrateoptics.raytracer.surface import Surface

from pyrateoptics.raytracer.aperture import CircularAperture
from pyrateoptics.raytracer.localcoordinates import LocalCoordinates
from pyrateoptics.raytracer.helpers import build_pilotbundle

from pyrateoptics import draw
from pyrateoptics.raytracer.globalconstants import degree

from pyrateoptics.raytracer.analysis.optical_system_analysis import\
    OpticalSystemAnalysis


logging.basicConfig(level=logging.DEBUG)

wavelength = 0.5876e-3

rnd_data1 = np.random.random((3, 3))  # np.eye(3)
rnd_data2 = np.random.random((3, 3))  # np.zeros((3, 3))#
lc = LocalCoordinates.p("1")
myeps = np.eye(3) + 0.1*rnd_data1 + 0.01*complex(0, 1)*rnd_data2
# aggressive complex choice of myeps
# myeps = np.eye(3) + 0.01*np.random.random((3, 3))
crystal = AnisotropicMaterial.p(lc, myeps)


# definition of optical system
s = OpticalSystem.p(matbackground=crystal)

lc0 = s.addLocalCoordinateSystem(
        LocalCoordinates.p(name="object", decz=0.0),
        refname=s.rootcoordinatesystem.name)
lc1 = s.addLocalCoordinateSystem(
        LocalCoordinates.p(name="m1", decz=50.0, tiltx=-math.pi/8),
        refname=lc0.name)  # objectDist
lc2 = s.addLocalCoordinateSystem(
        LocalCoordinates.p(name="m2_stop", decz=-50.0, decy=-20,
                         tiltx=math.pi/16), refname=lc1.name)
lc3 = s.addLocalCoordinateSystem(
        LocalCoordinates.p(name="m3", decz=50.0, decy=-30, tiltx=3*math.pi/32),
        refname=lc2.name)
lc4 = s.addLocalCoordinateSystem(
        LocalCoordinates.p(name="image1", decz=-50, decy=-15, tiltx=-math.pi/16),
        refname=lc3.name)
lc5 = s.addLocalCoordinateSystem(
        LocalCoordinates.p(name="oapara", decz=-100, decy=-35), refname=lc4.name)
lc5ap = s.addLocalCoordinateSystem(
        LocalCoordinates.p(name="oaparaap", decz=0, decy=35), refname=lc5.name)
lc6 = s.addLocalCoordinateSystem(
        LocalCoordinates.p(name="image2", decz=55, tiltx=1*math.pi/32),
        refname=lc5.name)

objectsurf = Surface.p(lc0)
m1surf = Surface.p(lc1, shape=Conic.p(lc1, curv=-0.01),
                 aperture=CircularAperture.p(lc1, maxradius=20.))
m2surf = Surface.p(lc2, shape=Conic.p(lc2, curv=0.01),
                 aperture=CircularAperture.p(lc2, maxradius=12.7))
m3surf = Surface.p(lc3, shape=Conic.p(lc3, curv=-0.006),
                 aperture=CircularAperture.p(lc3, maxradius=20.7))
image1 = Surface.p(lc4)
oapara = Surface.p(lc3, shape=Conic.p(lc5, curv=0.01, cc=-1.),
                 aperture=CircularAperture.p(lc5ap, maxradius=30.0))
image2 = Surface.p(lc6, aperture=CircularAperture.p(lc6, maxradius=20.0))


elem = OpticalElement.p(lc0, name="TMA")

# elem.addMaterial("crystal", crystal)

elem.addSurface("object", objectsurf, (None, None))
elem.addSurface("m1", m1surf, (None, None))
elem.addSurface("m2", m2surf, (None, None))
elem.addSurface("m3", m3surf, (None, None))
elem.addSurface("image1", image1, (None, None))
elem.addSurface("oapara", oapara, (None, None))
elem.addSurface("image2", image2, (None, None))

s.addElement("TMA", elem)

print(s.rootcoordinatesystem.pprint())

sysseq = [("TMA", [
            ("object", {}),
            ("m1", {"is_mirror": True}),
            ("m2", {"is_mirror": True}),
            ("m3", {"is_mirror": True}),
            ("image1", {}),
            ("oapara", {"is_mirror": True}),
            ("image2", {})])]

sysseq_pilot = [("TMA",
                 [
                    ("object", {}),
                    ("m1", {"is_mirror": True}),
                    ("m2", {"is_mirror": True}),
                    ("m3", {"is_mirror": True}),
                    ("m2", {"is_mirror": True}),
                    ("m1", {"is_mirror": True}),
                    ("m2", {"is_mirror": True}),
                    ("m1", {"is_mirror": True}),
                    ("m2", {"is_mirror": True})
                 ])
                ]


obj_dx = 0.1
obj_dphi = 1.*degree

osa = OpticalSystemAnalysis(s, sysseq, name="Analysis")

raysdict = {"opticalsystem": s, "startz": -5., "radius": 10.,
            "raster": raster.MeridionalFan()}
osa.aim(5, raysdict, bundletype="collimated", wave=wavelength)

r2 = osa.trace()[0]


kw = 5.*degree

pilotbundles = build_pilotbundle(objectsurf,
                                 crystal,
                                 (obj_dx, obj_dx),
                                 (obj_dphi, obj_dphi),
                                 kunitvector=np.array([0,
                                                       math.sin(kw),
                                                       math.cos(kw)]),
                                 num_sampling_points=3)
(pilotray2, r3) = s.para_seqtrace(pilotbundles[-1],
                                  osa.initial_bundles[0], sysseq)

draw(s, [(r2, "blue"), (r3, "orange"), (pilotray2, "red")])
