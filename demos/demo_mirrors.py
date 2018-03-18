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


from core import raster
from core.material_isotropic import ConstantIndexGlass
from core.material_anisotropic import AnisotropicMaterial
from core import surfShape
from core.optical_element import OpticalElement
from core.optical_element_analysis import OpticalElementAnalysis
from core.optical_system import OpticalSystem
from core.surface import Surface
from core.ray import RayBundle

from core.aperture import CircularAperture
from core.localcoordinates import LocalCoordinates

from core.globalconstants import canonical_ey

import math
import logging
logging.basicConfig(level=logging.DEBUG)

import core.helpers

wavelength = 0.5876e-3

# definition of optical system

#v = np.ones(3)# + 0.001*np.random.random(3)
#myeps = np.diag(v)


s = OpticalSystem() 

lc0 = s.addLocalCoordinateSystem(LocalCoordinates(name="object", decz=0.0), refname=s.rootcoordinatesystem.name)

#air = AnisotropicMaterial(lc0, myeps)  # tests for anisotropic mirror
air = ConstantIndexGlass(lc0, 1.0)
s.material_background = air

lc1 = s.addLocalCoordinateSystem(LocalCoordinates(name="m1", decz=50.0, tiltx=-math.pi/8), refname=lc0.name) # objectDist
lc2 = s.addLocalCoordinateSystem(LocalCoordinates(name="m2_stop", decz=-50.0, decy=-20, tiltx=math.pi/16), refname=lc1.name)
lc3 = s.addLocalCoordinateSystem(LocalCoordinates(name="m3", decz=50.0, decy=-30, tiltx=3*math.pi/32), refname=lc2.name)
lc4 = s.addLocalCoordinateSystem(LocalCoordinates(name="image1", decz=-50, decy=-15, tiltx=-math.pi/16), refname=lc3.name)
lc5 = s.addLocalCoordinateSystem(LocalCoordinates(name="oapara", decz=-100, decy=-35), refname=lc4.name)
lc5ap = s.addLocalCoordinateSystem(LocalCoordinates(name="oaparaap", decz=0, decy=35), refname=lc5.name)
lc6 = s.addLocalCoordinateSystem(LocalCoordinates(name="image2", decz=55, tiltx=1*math.pi/32), refname=lc5.name)
lc7 = s.addLocalCoordinateSystem(LocalCoordinates(name="image3", decz=5), refname=lc6.name)

objectsurf = Surface(lc0)
m1surf = Surface(lc1, shape=surfShape.Conic(lc1, curv=-0.01), apert=CircularAperture(lc1, 20.))
m2surf = Surface(lc2, shape=surfShape.Conic(lc2, curv=0.01), apert=CircularAperture(lc2, 12.7))
m3surf = Surface(lc3, shape=surfShape.Conic(lc3, curv=-0.006), apert=CircularAperture(lc3, 12.7))
image1 = Surface(lc4)
oapara = Surface(lc3, shape=surfShape.Conic(lc5, curv=0.01, cc=-1.), apert=CircularAperture(lc5ap, 30.0))
image2 = Surface(lc6, apert=CircularAperture(lc6, 20.0))
image3 = Surface(lc7, apert=CircularAperture(lc7, 20.0))


elem = OpticalElement(lc0, name="TMA")

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

rstobj = raster.MeridionalFan()
(px, py) = rstobj.getGrid(11)

rpup = 10
o = np.vstack((rpup*px, rpup*py, -5.*np.ones_like(px)))
k = np.zeros_like(o)
k[2,:] = 1.0 #2.*math.pi/wavelength
ey = np.zeros_like(o)
ey[1,:] =  1.
E0 = np.cross(k, ey, axisa=0, axisb=0).T

sysseq = [("TMA", 
           [
                ("object", {}), 
                ("m1", {"is_mirror":True}), 
                ("m2", {"is_stop":True, "is_mirror":True}), 
                ("m3", {"is_mirror":True}), 
                ("image1", {}), 
                ("oapara", {"is_mirror":True}), 
                ("image2", {}), 
                ("image3", {}) 
            ])
        ] 

# TODO: integrate mirrors (is_mirror:True) also into options dict

sysseq_pilot = [("TMA", 
                 [
                    ("object", True, {}), 
                    ("m1", False, {}), 
                    ("m2", False, {"is_stop":True}), 
                    ("m3", False, {}), 
                    ("m2", False, {}),
                    ("m1", False, {}),
                    ("m2", False, {}),
                    ("m1", False, {}),
                    ("m2", False, {})
                ])
                ] 
                
phi = 5.*math.pi/180.0

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


pilotbundles = core.helpers.build_pilotbundle(objectsurf, air, (obj_dx, obj_dx), (obj_dphi, obj_dphi), num_sampling_points=3)

rays_pilot = [s.seqtrace(p, sysseq) for p in pilotbundles]

(pilotray2, r3) = s.para_seqtrace(pilotbundles[-1], initialbundle, sysseq, use6x6=True)

(m_obj_stop, m_stop_img) = s.extractXYUV(pilotbundles[-1], sysseq, use6x6=True)

logging.info(np.array_str(m_obj_stop, precision=5, suppress_small=True))
logging.info(np.array_str(m_stop_img, precision=5, suppress_small=True))

logging.info(str(s.sequence_to_hitlist(sysseq)))


### TODO:
### first tries to implement aiming, but the code is somewhat hard to use
### we need to get rid of the pilot ray in every call
### we need to convert between XK representation local 3D coordinates and
### global raybundle coordinates in a more easy way
###
#oea = OpticalElementAnalysis(s.elements["TMA"])
#
#xyuvobjectstop = oea.calcXYUV([("object", "m1", 1), ("m1", "m2", 1)], pilotbundle2, sysseq[0][1], air)
#
#Axyuv = xyuvobjectstop[0:2, 0:2]
#Bxyuv = xyuvobjectstop[0:2, 2:4]
#Cxyuv = xyuvobjectstop[2:4, 0:2]
#Dxyuv = xyuvobjectstop[2:4, 2:4]
#
#
#
#alpha = np.linspace(0, 360, 20)*math.pi/180.
#pts = np.vstack((8.*np.cos(alpha), 8.*np.sin(alpha), np.zeros_like(alpha), np.zeros_like(alpha)))
#ptsXY = pts[0:2]
#
#kobj = np.dot(np.linalg.inv(Bxyuv), ptsXY)
#
#ptsobj = np.dot(np.linalg.inv(xyuvobjectstop), pts)
#
#
#
#xobj = ptsobj[0]
#yobj = ptsobj[1]
#
#kxobj = kobj[0]
#kyobj = kobj[1]
#
#o = np.vstack((xobj, yobj, np.zeros_like(xobj)))
#k = np.zeros_like(o)
#k[0,:] = kxobj
#k[1,:] = kyobj
#k[2,:] = np.sqrt((2.*math.pi/wavelength)**2 - kxobj**2 - kyobj**2)
#
#ey = np.zeros_like(o)
#ey[1,:] =  1.
#
#E0 = np.cross(k, ey, axisa=0, axisb=0).T
#initialbundle = RayBundle(x0=lc0.returnLocalToGlobalPoints(o), k0=lc0.returnLocalToGlobalDirections(k), Efield0=lc0.returnLocalToGlobalDirections(E0), wave=wavelength)
#r4 = s.seqtrace(initialbundle, sysseq)



fig = plt.figure(1)
ax = fig.add_subplot(111)
ax.axis('equal')
if StrictVersion(matplotlib.__version__) < StrictVersion('2.0.0'):
    ax.set_axis_bgcolor('white')
else:
    ax.set_facecolor('white')



phi = 0. #math.pi/4
pn = np.array([math.cos(phi), 0, math.sin(phi)]) # canonical_ex
up = canonical_ey

#print("drawing!")
for (i, r) in enumerate(r2):
    r.draw2d(ax, color="blue", plane_normal=pn, up=up)
for r_p in rays_pilot:
    for (i, r) in enumerate(r_p):    
        r.draw2d(ax, color="red", plane_normal=pn, up=up)

r3.draw2d(ax, color="orange", plane_normal=pn, up=up)
pilotray2.draw2d(ax, color="red", plane_normal=pn, up=up)

#r4.draw2d(ax, color="pink", plane_normal=pn, up=up)
#pilotray.draw2d(ax, color="darkgreen", plane_normal=pn, up=up)


s.draw2d(ax, color="grey", vertices=50, plane_normal=pn, up=up) # try for phi=0.
#s.draw2d(ax, color="grey", inyzplane=False, vertices=50, plane_normal=pn, up=up) # try for phi=pi/4

plt.show()

