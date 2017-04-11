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
import time
import math

from core import raster
from core.ray import RayBundle

from core.globalconstants import standard_wavelength

from core.helpers import build_simple_optical_system
from core.globalconstants import canonical_ey

wavelength = standard_wavelength

# definition of optical system

mat_dict = {"glass":1.7, "glass2":1.5}


(s, seq) = build_simple_optical_system(
        [(-5.922, 0, 2.0, "glass", "surf1"),
         (-3.160, 0, 3.0, None, "surf2"),
         (15.884, 0, 5.0, "glass", "surf3"),
        (-12.756, 0, 3.0, None, "surf4"),
        (0, 0, 3.0, None, "stop"),
        (3.125, 0, 2.0, "glass2", "surf5"),
        (1.479, 0, 3.0, None, "surf6"),
        (0, 0, 19.0, None, "surf7")
         ], mat_dict)

nrays = 100000
nrays_draw = 21


def collimated_bundle(nrays, startz, radius, rast):
    rstobj = rast
    (px, py) = rstobj.getGrid(nrays)
    rpup = radius
    o = np.vstack((rpup*px, rpup*py, startz*np.ones_like(px)))
    k = np.zeros_like(o)
    k[2,:] = 2.*math.pi/wavelength
    E0 = np.cross(k, canonical_ey, axisa=0, axisb=0).T
    return (o, k, E0)

# benchmark

# definition of rays
#nray = 1E5 # number of rays
#aimy = aim.aimFiniteByMakingASurfaceTheStop(s, pupilType=pupil.ObjectSpaceNA, #.StopDiameter,
#                                           pupilSizeParameter=0.2,#3.0,
#                                            fieldType= field.ObjectHeight,
#                                            rasterType= raster.RectGrid,
#                                            nray=nray, wavelength=wavelength, stopPosition=5)
#initialBundle = aimy.getInitialRayBundle(s, fieldXY=np.array([0, 0]), wavelength=wavelength)
#nray = len(initialBundle.o[0, :])

(x0, k0, E0) = collimated_bundle(nrays, -5., 1., raster.RectGrid())
t0 = time.clock()
initialraybundle = RayBundle(x0=x0, k0=k0, Efield0=E0)
raypath = s.seqtrace(initialraybundle, seq)
print "benchmark : ", time.clock() - t0, "s for tracing ", nrays, " rays through ", len(s.elements["stdelem"].surfaces) - 1, " surfaces."
print "             That is ", int(round(nrays * (len(s.elements["stdelem"].surfaces) - 1) / (time.clock() - t0))), "ray-surface-operations per second"

# plot

(x0_draw, k0_draw, E0_draw)= collimated_bundle(nrays_draw, -5., 1., raster.MeridionalFan())
initialraybundle_draw = RayBundle(x0=x0_draw, k0=k0_draw, Efield0=E0_draw)
raypath_draw = s.seqtrace(initialraybundle_draw, seq)


fig = plt.figure(1)
ax = fig.add_subplot(111)
ax.axis('equal')
ax.set_facecolor('white')


phi = 0.#math.pi/4
pn = np.array([math.cos(phi), 0, math.sin(phi)]) # canonical_ex
up = canonical_ey

raypath_draw.draw2d(ax, color="blue", plane_normal=pn, up=up) 
for e in s.elements.itervalues():
    for surfs in e.surfaces.itervalues():
        surfs.draw2d(ax, color="grey", vertices=50, plane_normal=pn, up=up) # try for phi=0.
        #surfs.draw2d(ax, color="grey", inyzplane=False, vertices=50, plane_normal=pn, up=up) # try for phi=pi/4

plt.show()


