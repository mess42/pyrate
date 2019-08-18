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

import math
import logging

import numpy as np

from pyrateoptics.sampling2d import raster
from pyrateoptics.raytracer.material.material_isotropic import\
    ConstantIndexGlass
from pyrateoptics.raytracer.surface_shape import Conic
from pyrateoptics.raytracer.optical_element import OpticalElement
from pyrateoptics.raytracer.optical_system import OpticalSystem
from pyrateoptics.raytracer.surface import Surface

from pyrateoptics.raytracer.aperture import CircularAperture
from pyrateoptics.raytracer.localcoordinates import LocalCoordinates

from pyrateoptics import draw

from pyrateoptics.raytracer.globalconstants import degree

from pyrateoptics.raytracer.aim import Aimy

logging.basicConfig(level=logging.INFO)


wavelength = 0.5876e-3

# definition of optical system

# v = np.ones(3)# + 0.001*np.random.random(3)
# myeps = np.diag(v)


s = OpticalSystem.p()

lc0 = s.addLocalCoordinateSystem(LocalCoordinates.p(name="object", decz=0.0),
                                 refname=s.rootcoordinatesystem.name)

# air = AnisotropicMaterial(lc0, myeps)  # tests for anisotropic mirror
air = ConstantIndexGlass.p(lc0, 1.0)
s.material_background = air

lc1 = s.addLocalCoordinateSystem(
            LocalCoordinates.p(name="m1", decz=50.0, tiltx=-math.pi/8),
            refname=lc0.name)  # objectDist
lc2 = s.addLocalCoordinateSystem(
            LocalCoordinates.p(name="m2_stop",
                               decz=-50.0,
                               decy=-20,
                               tiltx=math.pi/16),
            refname=lc1.name)
lc3 = s.addLocalCoordinateSystem(
            LocalCoordinates.p(name="m3", decz=50.0, decy=-30,
                               tiltx=3*math.pi/32), refname=lc2.name)
lc4 = s.addLocalCoordinateSystem(
            LocalCoordinates.p(name="image1", decz=-50, decy=-15,
                               tiltx=-math.pi/16), refname=lc3.name)
lc5 = s.addLocalCoordinateSystem(
            LocalCoordinates.p(name="oapara", decz=-100, decy=-35),
            refname=lc4.name)
lc5ap = s.addLocalCoordinateSystem(
            LocalCoordinates.p(name="oaparaap", decz=0, decy=35),
            refname=lc5.name)
lc6 = s.addLocalCoordinateSystem(
            LocalCoordinates.p(name="image2", decz=52.8, tiltx=1*math.pi/32),
            refname=lc5.name)
lc7 = s.addLocalCoordinateSystem(
            LocalCoordinates.p(name="image3", decz=5), refname=lc6.name)

objectsurf = Surface.p(lc0)
m1surf = Surface.p(lc1, shape=Conic.p(lc1, curv=-0.01))
m2surf = Surface.p(lc2, shape=Conic.p(lc2, curv=0.01))
m3surf = Surface.p(lc3, shape=Conic.p(lc3, curv=-0.006))
image1 = Surface.p(lc4)
oapara = Surface.p(lc3, shape=Conic.p(lc5, curv=0.01, cc=-1.))
image2 = Surface.p(lc6, aperture=CircularAperture(lc6, maxradius=20.0))
image3 = Surface.p(lc7, aperture=CircularAperture(lc7, maxradius=20.0))


elem = OpticalElement.p(lc0, name="TMA")

elem.addMaterial("air", air)

elem.addSurface("object", objectsurf, (None, None))
elem.addSurface("m1", m1surf, (None, None))
elem.addSurface("m2", m2surf, (None, None))
elem.addSurface("m3", m3surf, (None, None))
elem.addSurface("image1", image1, (None, None))
elem.addSurface("oapara", oapara, (None, None))
elem.addSurface("image2", image2, (None, None))
elem.addSurface("image3", image3, (None, None))

s.addElement("TMA", elem)

print(s.rootcoordinatesystem.pprint())

sysseq = [("TMA",
           [
                ("object", {}),
                ("m1", {"is_mirror": True}),
                ("m2", {"is_stop": True, "is_mirror": True}),
                ("m3", {"is_mirror": True}),
                ("image1", {}),
                ("oapara", {"is_mirror": True}),
                ("image2", {}),
                ("image3", {})
            ]
           )
          ]

a = Aimy(s, sysseq, name="Aimy", stopsize=2., num_pupil_points=5)
a.pupil_raster = raster.MeridionalFan()


def correctKRayBundle(bundle):
    """
    Should get correct k from raybundle.
    """
    pass


initbundle1 = a.aim(np.array([0, 0]))
initbundle2 = a.aim(np.array([0, 0.5*degree]))
initbundle3 = a.aim(np.array([0, -0.5*degree]))

(pp1, r1p) = s.para_seqtrace(a.pilotbundle, initbundle1, sysseq)
(pp2, r2p) = s.para_seqtrace(a.pilotbundle, initbundle2, sysseq)
(pp3, r3p) = s.para_seqtrace(a.pilotbundle, initbundle3, sysseq)

r1r = s.seqtrace(initbundle1, sysseq)
r2r = s.seqtrace(initbundle2, sysseq)
r3r = s.seqtrace(initbundle3, sysseq)


draw(s, [(r1p, "blue"), (r2p, "green"), (r3p, "orange")])
# draw(s, [(r1r, "blue"), (r2r, "green"), (r3r, "orange")])


# TODO:
# first tries to implement aiming, but the code is somewhat hard to use
# we need to get rid of the pilot ray in every call
# we need to convert between XK representation local 3D coordinates and
# global raybundle coordinates in a more easy way
#
# oea = OpticalElementAnalysis(s.elements["TMA"])
#
# xyuvobjectstop = oea.calcXYUV([("object", "m1", 1), ("m1", "m2", 1)],
#                               pilotbundle2, sysseq[0][1], air)
#
# Axyuv = xyuvobjectstop[0:2, 0:2]
# Bxyuv = xyuvobjectstop[0:2, 2:4]
# Cxyuv = xyuvobjectstop[2:4, 0:2]
# Dxyuv = xyuvobjectstop[2:4, 2:4]
#
#
#
# alpha = np.linspace(0, 360, 20)*math.pi/180.
# pts = np.vstack((8.*np.cos(alpha), 8.*np.sin(alpha), np.zeros_like(alpha),
#                 np.zeros_like(alpha)))
# ptsXY = pts[0:2]
#
# kobj = np.dot(np.linalg.inv(Bxyuv), ptsXY)
#
# ptsobj = np.dot(np.linalg.inv(xyuvobjectstop), pts)
#
#
#
# xobj = ptsobj[0]
# yobj = ptsobj[1]
#
# kxobj = kobj[0]
# kyobj = kobj[1]
#
# o = np.vstack((xobj, yobj, np.zeros_like(xobj)))
# k = np.zeros_like(o)
# k[0,:] = kxobj
# k[1,:] = kyobj
# k[2,:] = np.sqrt((2.*math.pi/wavelength)**2 - kxobj**2 - kyobj**2)
#
# ey = np.zeros_like(o)
# ey[1,:] =  1.
#
# E0 = np.cross(k, ey, axisa=0, axisb=0).T
# initialbundle = RayBundle(x0=lc0.returnLocalToGlobalPoints(o),
#                           k0=lc0.returnLocalToGlobalDirections(k),
#                           Efield0=lc0.returnLocalToGlobalDirections(E0),
#                           wave=wavelength)
# r4 = s.seqtrace(initialbundle, sysseq)
