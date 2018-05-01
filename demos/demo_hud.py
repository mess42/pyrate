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
import logging


from pyrateoptics.sampling2d import raster
from pyrateoptics.material.material_isotropic import ConstantIndexGlass, ModelGlass
from pyrateoptics.material.material_anisotropic import AnisotropicMaterial
from pyrateoptics.raytracer.surfShape import Conic, Biconic
from pyrateoptics.raytracer.optical_element import OpticalElement
from pyrateoptics.analysis.optical_element_analysis import OpticalElementAnalysis
from pyrateoptics.raytracer.optical_system import OpticalSystem
from pyrateoptics.raytracer.surface import Surface
from pyrateoptics.raytracer.ray import RayBundle

from pyrateoptics.raytracer.aperture import CircularAperture
from pyrateoptics.raytracer.localcoordinates import LocalCoordinates

from pyrateoptics.raytracer.globalconstants import canonical_ey, degree, standard_wavelength

import math

import pyrateoptics.raytracer.helpers

# definition of optical system

# Design: US patent no. 5701202 A, inventor: Koichi Takahashi
# and also: Bo Chen, Alois M. Herkommer, Opt. Express 24, 26999 (2016)

logging.basicConfig(level=logging.INFO)

s = OpticalSystem() 

lc0 = s.addLocalCoordinateSystem(LocalCoordinates(name="object", decz=0.0), refname=s.rootcoordinatesystem.name)

air = ConstantIndexGlass(lc0, 1.0)
glass = ConstantIndexGlass(lc0, 1.492)
s.material_background = air

si = -1.

lcD1 = s.addLocalCoordinateSystem(LocalCoordinates(name="D1", decz=30.002), refname=lc0.name)
lcS1 = s.addLocalCoordinateSystem(LocalCoordinates(name="S1",               decy=-24.028,   decz=26.360, tiltx=-si*14.7*degree, tiltThenDecenter=False), refname=lc0.name)
lcD1prime = s.addLocalCoordinateSystem(LocalCoordinates(name="D1prime",     decy=0.,        decz=30.002, tiltx=-si*1.066*degree, tiltThenDecenter=False), refname=lc0.name)
lcD2 = s.addLocalCoordinateSystem(LocalCoordinates(name="D2",               decy=-0.251,    decz=43.485, tiltx=-si*1.066*degree, tiltThenDecenter=False), refname=lc0.name)
lcS2 = s.addLocalCoordinateSystem(LocalCoordinates(name="S2",               decy=19.109,    decz=33.339, tiltx=si*36.660*degree, tiltThenDecenter=False), refname=lc0.name)
lcD2prime = s.addLocalCoordinateSystem(LocalCoordinates(name="D2prime",     decy=-0.251,    decz=43.485, tiltx=si*38.376*degree, tiltThenDecenter=False), refname=lc0.name)
lcD3 = s.addLocalCoordinateSystem(LocalCoordinates(name="D3",               decy=-11.858,   decz=28.827, tiltx=si*38.376*degree, tiltThenDecenter=False), refname=lc0.name)
lcS3 = s.addLocalCoordinateSystem(LocalCoordinates(name="S3",               decy=-24.028,   decz=26.360, tiltx=-si*14.7*degree, tiltThenDecenter=False), refname=lc0.name)
lcD3prime = s.addLocalCoordinateSystem(LocalCoordinates(name="D3prime",     decy=-11.858,   decz=28.827, tiltx=-si*55.019*degree, tiltThenDecenter=False), refname=lc0.name)
lcD4 = s.addLocalCoordinateSystem(LocalCoordinates(name="D4",               decy=-23.067,   decz=36.667, tiltx=-si*55.019*degree, tiltThenDecenter=False), refname=lc0.name)
lcS4 = s.addLocalCoordinateSystem(LocalCoordinates(name="S4",               decy=-35.215,   decz=18.817, tiltx=-si*47.770*degree, tiltThenDecenter=False), refname=lc0.name)
lcD4prime = s.addLocalCoordinateSystem(LocalCoordinates(name="D4prime",     decy=-23.067,   decz=36.667, tiltx=-si*50.668*degree, tiltThenDecenter=False), refname=lc0.name)
lcimage = s.addLocalCoordinateSystem(LocalCoordinates(name="image",         decy=-30.892,   decz=43.083, tiltx=-si*50.668*degree, tiltThenDecenter=False), refname=lc0.name)

objsurf = Surface(lc0)
D1surf = Surface(lcD1)
#S1surf = Surface(lcS1, shape=Conic(lcS1, curv=si*1./108.187, cc=0), apert=CircularAperture(lcS1, 100.0))
S1surf = Surface(lcS1, shape=Biconic(lcS1, curvy=si*1./108.187, curvx=si*1./73.105, coefficients=[(0., 0.), (-si*5.542e-7, -0.08), (-si*8.176e-11, -1.379)]), apert=CircularAperture(lcS1, 40.0))
D1Psurf = Surface(lcD1prime)
D2surf = Surface(lcD2)
#S2surf = Surface(lcS2, shape=Conic(lcS2, curv=si*1./69.871, cc=-0.1368), apert=CircularAperture(lcS2, 60.0))
S2surf = Surface(lcS2, shape=Biconic(lcS2, curvy=si*1./69.871, curvx=si*1./60.374, ccy=-0.1368, ccx=-0.123, coefficients=[(0., 0.), (si*7.233e-11, 29.075), (si*4.529e-12, -2.085)]), apert=CircularAperture(lcS2, 40.0))
D2Psurf = Surface(lcD2prime)
D3surf = Surface(lcD3)
#S3surf = Surface(lcS3, shape=Conic(lcS3, curv=si*1./108.187, cc=0.0001), apert=CircularAperture(lcS3, 100.0))
S3surf = Surface(lcS3, shape=Biconic(lcS3, curvy=si*1./108.187, curvx=si*1./73.105, coefficients=[(0., 0.), (-si*5.542e-7, -0.08), (-si*8.176e-11, -1.379)]), apert=CircularAperture(lcS3, 40.0))
D3Psurf = Surface(lcD3prime)
D4surf = Surface(lcD4)
S4surf = Surface(lcS4, shape=Conic(lcS4, curv=1./77.772), apert=CircularAperture(lcS4, 40.0))
D4Psurf = Surface(lcD4prime)
imgsurf = Surface(lcimage)

elem = OpticalElement(lc0, name="HUD")

elem.addMaterial("air", air)
elem.addMaterial("glass", glass)

elem.addSurface("object", objsurf, (None, None))
elem.addSurface("d1", D1surf, (None, None))
elem.addSurface("s1", S1surf, (None, "glass"))
elem.addSurface("d1p", D1Psurf, ("glass", "glass"))
elem.addSurface("d2", D2surf, ("glass", "glass"))
elem.addSurface("s2", S2surf, ("glass", "glass"))
elem.addSurface("d2p", D2Psurf, ("glass", "glass"))
elem.addSurface("d3", D3surf, ("glass", "glass"))
elem.addSurface("s3", S1surf, ("glass", "glass"))
elem.addSurface("d3p", D3Psurf, ("glass", "glass"))
elem.addSurface("d4", D4surf, ("glass", "glass"))
elem.addSurface("s4", S4surf, ("glass", None))
elem.addSurface("d4p", D4Psurf, (None, None))
elem.addSurface("image", imgsurf, (None, None))

s.addElement("HUD", elem)

print(s.rootcoordinatesystem.pprint())

rstobj = raster.MeridionalFan()
(px, py) = rstobj.getGrid(3)

rpup = 2
o = np.vstack((rpup*px, rpup*py, 0.*np.ones_like(px)))
k1 = np.zeros_like(o)
k1[2,:] = 1.0 
ey = np.zeros_like(o)
ey[1,:] =  1.
E1 = np.cross(k1, ey, axisa=0, axisb=0).T

k2 = np.zeros_like(o)
k2[1,:] = math.sin(15*degree) 
k2[2,:] = math.cos(15*degree) 
E2 = np.cross(k2, ey, axisa=0, axisb=0).T

k3 = np.zeros_like(o)
k3[1,:] = math.sin(-15*degree) 
k3[2,:] = math.cos(-15*degree)
E3 = np.cross(k3, ey, axisa=0, axisb=0).T


sysseq = [("HUD", 
           [
                ("object", {"is_stop":True}), 
                ("d1", {}), 
                ("s1", {}), 
                ("d1p", {}), 
                ("d2", {}), 
                ("s2", {"is_mirror":True}), 
                ("d2p", {}), 
                ("d3", {}), 
                ("s3", {"is_mirror":True}), 
                ("d3p", {}), 
                ("d4", {}), 
                ("s4", {}), 
                ("d4p", {}), 
                ("image", {}) 
            ])
        ] 


initialbundle1 = RayBundle(x0=o, k0=k1, Efield0=E1, wave=standard_wavelength)
initialbundle2 = RayBundle(x0=o, k0=k2, Efield0=E2, wave=standard_wavelength)
initialbundle3 = RayBundle(x0=o, k0=k3, Efield0=E3, wave=standard_wavelength)
r1 = s.seqtrace(initialbundle1, sysseq)
r2 = s.seqtrace(initialbundle2, sysseq)
r3 = s.seqtrace(initialbundle3, sysseq)

obj_dx = 0.1
obj_dphi = 5*degree
pilotbundles = pyrateoptics.raytracer.helpers.build_pilotbundle(objsurf, air, (obj_dx, obj_dx), (obj_dphi, obj_dphi), num_sampling_points=3)

rays_pilot = [s.seqtrace(p, sysseq) for p in pilotbundles[2:]]
# only last two bundles hit the next surface


(pilotray, r_pilot) = s.para_seqtrace(pilotbundles[-1], initialbundle1, sysseq, use6x6=True)

(m_obj_stop, m_stop_img) = s.extractXYUV(pilotbundles[-1], sysseq, use6x6=True)

print(np.array_str(m_obj_stop, precision=5, suppress_small=True))
print(np.array_str(m_stop_img, precision=5, suppress_small=True))

#print(s.sequence_to_hitlist(sysseq))




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
for (i, r) in enumerate(r1):
    r.draw2d(ax, color="red", plane_normal=pn, up=up)
for (i, r) in enumerate(r3):
    r.draw2d(ax, color="green", plane_normal=pn, up=up)
for r_p in rays_pilot:
    for (i, r) in enumerate(r_p):    
        r.draw2d(ax, color="orange", plane_normal=pn, up=up)

pz = np.array([26.36, 33.339, 18.817])
py = np.array([-24.028, 19.109, -35.215])

ax.scatter(pz, py)

s.draw2d(ax, color="grey", vertices=50, plane_normal=pn, up=up) # try for phi=0.
#s.draw2d(ax, color="grey", inyzplane=False, vertices=50, plane_normal=pn, up=up) # try for phi=pi/4

plt.show()

