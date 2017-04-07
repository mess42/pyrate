#!/usr/bin/env/python
"""
Pyrate - Optical raytracing based on Python

Copyright (C) 2017 Moritz Esslinger <moritz.esslinger@web.de>
               and Johannes Hartung <j.hartung@gmx.net>
               and     Uwe Lippmann <uwe.lippmann@web.de>
               and    Thomas Heinze <t.heinze@fn.de>

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

from core import raster
from core import material
from core import surfShape
from core.optical_element import OpticalElement
from core.optical_element_analysis import OpticalElementAnalysis
from core.optical_system import OpticalSystem
from core.surface import Surface
from core.ray import RayBundle

from core.aperture import CircularAperture
from core.localcoordinates import LocalCoordinates

from core.globalconstants import canonical_ey

from core.helpers import sequence_to_hitlist

import math

wavelength = 0.5876e-3

# definition of optical system
s = OpticalSystem() 

lc0 = s.addLocalCoordinateSystem(LocalCoordinates(name="object", decz=0.0), refname=s.rootcoordinatesystem.name)
lc1 = s.addLocalCoordinateSystem(LocalCoordinates(name="m1", decz=50.0, tiltx=-math.pi/8), refname=lc0.name) # objectDist
lc2 = s.addLocalCoordinateSystem(LocalCoordinates(name="m2_stop", decz=-50.0, decy=-20, tiltx=math.pi/16), refname=lc1.name)
lc3 = s.addLocalCoordinateSystem(LocalCoordinates(name="m3", decz=50.0, decy=-30, tiltx=math.pi/8), refname=lc2.name)
lc4 = s.addLocalCoordinateSystem(LocalCoordinates(name="image1", decz=-50, decy=-15, tiltx=-math.pi/16), refname=lc3.name)
lc5 = s.addLocalCoordinateSystem(LocalCoordinates(name="oapara", decz=-100, decy=-35), refname=lc4.name)
lc5ap = s.addLocalCoordinateSystem(LocalCoordinates(name="oaparaap", decz=0, decy=35), refname=lc5.name)
lc6 = s.addLocalCoordinateSystem(LocalCoordinates(name="image2", decz=55), refname=lc5.name)

stopsurf = Surface(lc0)
frontsurf = Surface(lc1, shape=surfShape.Conic(lc1, curv=-0.001), apert=CircularAperture(lc1, 20.))
cementsurf = Surface(lc2, shape=surfShape.Conic(lc2, curv=-0.001), apert=CircularAperture(lc2, 12.7))
rearsurf = Surface(lc3, shape=surfShape.Conic(lc3, curv=0.001), apert=CircularAperture(lc3, 12.7))
image1 = Surface(lc4)
oapara = Surface(lc3, shape=surfShape.Conic(lc5, curv=0.01, cc=-1.), apert=CircularAperture(lc5ap, 12.7))
image2 = Surface(lc6)

air = material.ConstantIndexGlass(lc0, 1.0)

elem = OpticalElement(lc0, label="TMA")

elem.addMaterial("air", air)

elem.addSurface("object", stopsurf, (None, None))
elem.addSurface("m1", frontsurf, (None, None))
elem.addSurface("m2", cementsurf, (None, None))
elem.addSurface("m3", rearsurf, (None, None))
elem.addSurface("image1", image1, (None, None))
elem.addSurface("oapara", oapara, (None, None))
elem.addSurface("image2", image2, (None, None))

s.addElement("TMA", elem)

print(s.rootcoordinatesystem.pprint())

rstobj = raster.MeridionalFan()
(px, py) = rstobj.getGrid(11)

rpup = 10
o = np.vstack((rpup*px, rpup*py, -5.*np.ones_like(px)))

k = np.zeros_like(o)
k[2,:] = 2.*math.pi/wavelength

ey = np.zeros_like(o)
ey[1,:] =  1.

E0 = np.cross(k, ey, axisa=0, axisb=0).T

sysseq = [("TMA", [("object", True, True), ("m1", False, True), ("m2", False, True), ("m3", False, True), ("image1", True, True), ("oapara", False, True), ("image2", True, True) ])] 

sysseq_pilot = [("TMA", 
                 [
                    ("object", True, True), 
                    ("m1", False, True), 
                    ("m2", False, True), 
                    ("m3", False, True), 
                    ("m2", False, True),
                    ("m1", False, True),
                    ("m2", False, True),
                    ("m1", False, True),
                    ("m2", False, True)
                ])
                ] 
                
print(sequence_to_hitlist(sysseq_pilot[0][1]))


phi = 5.*math.pi/180.0

obj_dx = 0.1
obj_dphi = 1.*math.pi/180.0

kwave = 2.*math.pi/wavelength

initialbundle = RayBundle(x0=o, k0=k, Efield0=E0, wave=wavelength)
r2 = s.seqtrace(initialbundle, sysseq)

pilotbundle = RayBundle(
                x0 = np.array([[0], [0], [0]]), 
                k0 = np.array([[0], [kwave*math.sin(phi)], [kwave*math.cos(phi)]]), 
                Efield0 = np.array([[1], [0], [0]]), wave=wavelength
                )

pilotbundle2 = RayBundle(
                x0 = np.array([[0, obj_dx, 0, 0, 0], [0, 0, obj_dx, 0, 0], [0, 0, 0, 0, 0]]), 
                k0 = np.array([[0, 0, 0, kwave*math.sin(obj_dphi), 0], [0, 0, 0, 0, kwave*math.sin(obj_dphi)], [kwave, kwave, kwave, kwave*math.cos(obj_dphi), kwave*math.cos(obj_dphi)]]), 
                Efield0 = np.array([[1, 1, 1, 1, 1], [0, 0, 0, 0, 0], [0, 0, 0, 0, 0]]), wave=wavelength
                )


pilotray = s.seqtrace(pilotbundle, sysseq_pilot)
#pilotray2 = s.elements["TMA"].seqtrace(pilotbundle2, sysseq[0][1], air)

oea = OpticalElementAnalysis(s.elements["TMA"])

(pilotray2, matrices) = oea.calculateXYUV(pilotbundle2, sysseq[0][1], air)

# TODO: these dictionary matrices are not compatible with bouncing pilotray and his buddies
# TODO: tracing from object to image2 (where the focus is), picture looks good

mat = np.dot(matrices[('image2', 'oapara')], 
             np.dot(matrices[('oapara', 'image1')], 
                    np.dot(matrices[('image1', 'm3')], 
                           np.dot(matrices[('m3', 'm2')], 
                                  np.dot(matrices[('m2', 'm1')], 
                                         matrices[('m1', 'object')])))))
alpha = np.arange(0, 360, 10)

pts = np.dot(mat, np.vstack((np.cos(alpha*math.pi/180.0), np.sin(alpha*math.pi/180.0), np.zeros_like(alpha), np.zeros_like(alpha))))
xf = pts[0]
yf = pts[1]

fig = plt.figure(1)
ax = fig.add_subplot(121)
ax2 = fig.add_subplot(122)

ax.axis('equal')
ax2.axis('equal')

ax.set_axis_bgcolor('white')



phi = 0. #math.pi/4
pn = np.array([math.cos(phi), 0, math.sin(phi)]) # canonical_ex
up = canonical_ey

ax2.plot(xf, yf, "ro")

print("drawing!")
r2.draw2d(ax, color="blue", plane_normal=pn, up=up) 
pilotray.draw2d(ax, color="darkgreen", plane_normal=pn, up=up)
pilotray2.draw2d(ax, color="red", plane_normal=pn, up=up)
for e in s.elements.itervalues():
    for surfs in e.surfaces.itervalues():
        surfs.draw2d(ax, color="grey", vertices=50, plane_normal=pn, up=up) # try for phi=0.
        #surfs.draw2d(ax, color="grey", inyzplane=False, vertices=50, plane_normal=pn, up=up) # try for phi=pi/4

plt.show()

