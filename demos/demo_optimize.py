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

import sys
import logging
import math
from distutils.version import StrictVersion


import numpy as np
import matplotlib.pyplot as plt
import matplotlib

from pyrateoptics import listOptimizableVariables
from pyrateoptics.raytracer.material.material_isotropic import\
    ConstantIndexGlass
from pyrateoptics.raytracer.surface_shape import Conic
from pyrateoptics.raytracer.ray import RayBundle

from pyrateoptics.raytracer.localcoordinates import LocalCoordinates

from pyrateoptics.raytracer.globalconstants import standard_wavelength

from pyrateoptics.raytracer.optical_element import OpticalElement
from pyrateoptics.raytracer.optical_system import OpticalSystem
from pyrateoptics.raytracer.surface import Surface

from pyrateoptics.raytracer.globalconstants import canonical_ey, degree

from pyrateoptics.raytracer.analysis.optical_system_analysis import\
    OpticalSystemAnalysis
from pyrateoptics.raytracer.analysis.surface_shape_analysis import\
    ShapeAnalysis

from pyrateoptics.sampling2d.raster import RandomGrid

from pyrateoptics.optimize.optimize import Optimizer
from pyrateoptics.optimize.optimize_backends import (ScipyBackend,
                                                     Newton1DBackend,
                                                     ParticleSwarmBackend,
                                                     SimulatedAnnealingBackend)

wavelength = standard_wavelength

logging.basicConfig(level=logging.DEBUG)

# definition of optical system
s = OpticalSystem.p()  # objectDistance = 2.0

lc0 = s.addLocalCoordinateSystem(
            LocalCoordinates.p(name="object", decz=0.0),
            refname=s.rootcoordinatesystem.name)
lc1 = s.addLocalCoordinateSystem(
            LocalCoordinates.p(name="surf1", decz=2.0), refname=lc0.name)
# objectDist
lc2 = s.addLocalCoordinateSystem(
            LocalCoordinates.p(name="surf2", decz=3.0), refname=lc1.name)
lc3 = s.addLocalCoordinateSystem(
            LocalCoordinates.p(name="surf3", decz=5.0, tiltx=2.5*degree),
            refname=lc2.name)
lc4 = s.addLocalCoordinateSystem(
            LocalCoordinates.p(name="surf4", decz=3.0), refname=lc3.name)
lc5 = s.addLocalCoordinateSystem(
            LocalCoordinates.p(name="surf5", decz=3.0), refname=lc4.name)
lc6 = s.addLocalCoordinateSystem(
            LocalCoordinates.p(name="surf6", decz=2.0), refname=lc5.name)
lc7 = s.addLocalCoordinateSystem(
            LocalCoordinates.p(name="surf7", decz=3.0), refname=lc6.name)
lc8 = s.addLocalCoordinateSystem(
            LocalCoordinates.p(name="image", decz=19.0), refname=lc7.name)

objectsurf = Surface.p(lc0)
surf1 = Surface.p(lc1, shape=Conic.p(lc1, curv=1/-5.922))
surf2 = Surface.p(lc2, shape=Conic.p(lc2, curv=1/-3.160))
surf3 = Surface.p(lc3, shape=Conic.p(lc3, curv=1/15.884))
surf4 = Surface.p(lc4, shape=Conic.p(lc4, curv=1/-12.756))
stopsurf = Surface.p(lc5)
surf6 = Surface.p(lc6, shape=Conic.p(lc6, curv=1/3.125))
surf7 = Surface.p(lc7, shape=Conic.p(lc7, curv=0.1*1/1.479))
image = Surface.p(lc8)


elem = OpticalElement.p(lc0, name="lenssystem")

glass = ConstantIndexGlass.p(lc0, n=1.7)
glass2 = ConstantIndexGlass.p(lc0, n=1.5)

elem.addMaterial("glass", glass)
elem.addMaterial("glass2", glass2)

elem.addSurface("object", objectsurf, (None, None))
elem.addSurface("surf1", surf1, (None, "glass"))
elem.addSurface("surf2", surf2, ("glass", None))
elem.addSurface("surf3", surf3, (None, "glass"))
elem.addSurface("surf4", surf4, ("glass", None))
elem.addSurface("stop", stopsurf, (None, None))
elem.addSurface("surf6", surf6, (None, "glass2"))
elem.addSurface("surf7", surf7, ("glass2", None))
elem.addSurface("image", image, (None, None))

s.addElement("lenssys", elem)


# plot initial system
# aimy = aim.aimFiniteByMakingASurfaceTheStop(s, pupilType=pupil.ObjectSpaceNA,
#                                            pupilSizeParameter=0.2,#3.0,
#                                            fieldType= field.ObjectHeight,
#                                            rasterType= raster.RectGrid,
#                                            nray=20, wavelength=wavelength,
#                                            stopPosition=5)
#
# initialBundle2 = aimy.getInitialRayBundle(s, fieldXY=np.array([0, 0]),
#                                           wavelength=wavelength)
#
# r2 = RayPath(initialBundle2, s)
#
# initialBundle3 = aimy.getInitialRayBundle(s, fieldXY=np.array([0, 0.1]),
#                                           wavelength=wavelength)
# r3 = RayPath(initialBundle3, s)
#
# initialBundle4 = aimy.getInitialRayBundle(s, fieldXY=np.array([0, -0.1]),
#                                           wavelength=wavelength)
# r4 = RayPath(initialBundle4, s)
#
# fig = plt.figure(1)
# ax = fig.add_subplot(211)
# ax2 = fig.add_subplot(212)
#
# ax.axis('equal')
# ax.set_facecolor('black')
# ax2.axis('equal')
# ax2.set_facecolor('black')
#
# plots.drawLayout2d(ax, s, [r2, r3, r4])


# optimize
# print("Initial   merit function: ",
#       merit.mySimpleDumbRMSSpotSizeMeritFunction(s))


# reintroduced apertures after optimization run
# s.surfaces[1].aperture = CircularAperture(maxradius=0.55)
# s.surfaces[2].aperture = CircularAperture(maxradius=1.0)
# s.surfaces[3].aperture = CircularAperture(maxradius=1.3)
# s.surfaces[4].aperture = CircularAperture(maxradius=1.3)
# s.surfaces[5].aperture = CircularAperture(maxradius=1.01)
# s.surfaces[6].aperture = CircularAperture(maxradius=1.0)
# s.surfaces[7].aperture = CircularAperture(maxradius=1.0)


# aimy.setPupilRaster(rasterType= raster.RectGrid, nray=100)
# initialBundle2 = aimy.getInitialRayBundle(s, fieldXY=np.array([0, 0]),
#                                           wavelength=wavelength)
# r2 = RayPath(initialBundle2, s)
# initialBundle3 = aimy.getInitialRayBundle(s, fieldXY=np.array([0, 0.1]),
#                                           wavelength=wavelength)
# r3 = RayPath(initialBundle3, s)
# initialBundle4 = aimy.getInitialRayBundle(s, fieldXY=np.array([0, -0.1]),
#                                           wavelength=wavelength)
# r4 = RayPath(initialBundle4, s)


# fig15 = plt.figure(15)
# ax3 = fig15.add_subplot(111)
#
# plt.subplots_adjust(left=0.0, right=1.0, bottom=0.0, top=1.0)
# plots.drawLayout2d(ax2, s, [r2, r3, r4])
# plt.subplots_adjust(left=0.0, right=1.0, bottom=0.0, top=1.0)
# plots.drawSpotDiagram(ax3, s, r3, -1)
# plt.subplots_adjust(left=0.0, right=1.0, bottom=0.0, top=1.0)


sysseq = [("lenssys", [
            ("object", {}),
            ("surf1", {}),
            ("surf2", {}),
            ("surf3", {}),
            ("surf4", {}),
            ("stop", {"is_stop": True}),
            ("surf6", {}),
            ("surf7", {}),
            ("image", {})])]


osa = OpticalSystemAnalysis(s, sysseq, name="Analysis")

divbundledict = {"radius": 10*degree, "raster": RandomGrid()}
(o, k, E0) = osa.divergent_bundle(900, divbundledict, wave=wavelength)
initialbundle = RayBundle(x0=o, k0=k, Efield0=E0)
r2 = s.seqtrace(initialbundle, sysseq)

fig = plt.figure(1)
ax = fig.add_subplot(311)
ax2 = fig.add_subplot(312)
ax3 = fig.add_subplot(313)

ax.axis('equal')
ax2.axis('equal')
if StrictVersion(matplotlib.__version__) < StrictVersion('2.0.0'):
    ax.set_axis_bgcolor('white')
else:
    ax.set_facecolor('white')

phi = 0.  # math.pi/4
pn = np.array([math.cos(phi), 0, math.sin(phi)])  # canonical_ex
up = canonical_ey

for r in r2:
    r.draw2d(ax, color="blue", plane_normal=pn, up=up)
s.draw2d(ax, color="grey", vertices=50, plane_normal=pn, up=up)  # try phi=0.
# s.draw2d(ax, color="grey", inyzplane=False, vertices=50, plane_normal=pn,
#          up=up) # try for phi=pi/4

curv2 = s.elements["lenssys"].surfaces["surf2"].shape.curvature
curv2.toVariable()
curv2.setInterval(left=-0.35, right=0.35)
curv3 = s.elements["lenssys"].surfaces["surf3"].shape.curvature
curv3.toVariable()
curv3.setInterval(left=-0.35, right=0.35)
curv4 = s.elements["lenssys"].surfaces["surf4"].shape.curvature
curv4.toVariable()
curv4.setInterval(left=-0.35, right=0.35)
curv6 = s.elements["lenssys"].surfaces["surf6"].shape.curvature
curv6.toVariable()
curv6.setInterval(left=-0.35, right=0.35)
tltx_var = s.elements["lenssys"].surfaces["surf3"].rootcoordinatesystem.tiltx
tltx_var.toVariable()
tltx_var.setInterval(left=-3.*math.pi/180., right=3.*math.pi/180.)

listOptimizableVariables(s, filter_status='variable', max_line_width=80)


def osnone(my_s):
    """
    Do nothing
    """
    pass


def osupdate(my_s):
    """
    Update all coordinate systems during run
    """
    my_s.rootcoordinatesystem.update()


def meritfunctionrms(my_s):
    """
    Merit function for tracing a raybundle through system and calculate
    rms spot radius without centroid subtraction. Punish low number of
    rays, too.
    """
    my_initialbundle = RayBundle(x0=o, k0=k, Efield0=E0, wave=wavelength)
    rpaths = my_s.seqtrace(my_initialbundle, sysseq)

    x = rpaths[0].raybundles[-1].x[-1, 0, :]
    y = rpaths[0].raybundles[-1].x[-1, 1, :]

    xmean = np.mean(x)
    ymean = np.mean(y)

    res = np.sum((x - xmean)**2 + (y - ymean)**2) + 10.*math.exp(-len(x))

    return res


opt_backend = ScipyBackend(method='Nelder-Mead',
                           options={'maxiter': 1000, 'disp': True}, tol=1e-8)
if len(sys.argv) > 1:
    first_arg = sys.argv[1]

    if first_arg == "--help":
        print("p ... Particle Swarm (c1, c2)")
        print("s ... Simulated Annealing (Nt, Tt)")
        print("n ... Newtonian 1D (dx, iterations)")
        print("m ... Nelder-Mead")
        sys.exit(0)
    elif first_arg == "p":
        c1 = 2.2
        c2 = 2.1
        if len(sys.argv) == 4:
            (c1, c2) = (float(sys.argv[2]), float(sys.argv[3]))
        opt_backend = ParticleSwarmBackend(c1=2.2, c2=2.1)
    elif first_arg == "s":
        Nt = 30
        T0 = 20.
        if len(sys.argv) == 4:
            (Nt, T0) = (int(sys.argv[2]), float(sys.argv[3]))
        opt_backend = SimulatedAnnealingBackend(Nt=Nt,
                                                Tt=T0*np.exp(-np.linspace(0, 10, 20)))
    elif first_arg == "n":
        dx = 1e-6
        iterations = 100
        if len(sys.argv) == 4:
            (dx, iterations) = (float(sys.argv[2]), int(sys.argv[3]))
        opt_backend = Newton1DBackend(dx=dx, iterations=iterations)

# pylint3 E1130: false positive

optimi = Optimizer(s,
                   meritfunctionrms,
                   backend=opt_backend,
                   updatefunction=osupdate)
optimi.logger.setLevel(logging.DEBUG)
s = optimi.run()

r2 = s.seqtrace(initialbundle, sysseq)  # trace again

for r in r2:
    r.draw2d(ax2, color="blue", plane_normal=pn, up=up)

listOptimizableVariables(s, filter_status='variable', max_line_width=80)

s.draw2d(ax2, color="grey", vertices=50, plane_normal=pn, up=up)  # try phi=0.
# s.draw2d(ax, color="grey", inyzplane=False, vertices=50, plane_normal=pn, up=up) # try for phi=pi/4

osa.aim(500, divbundledict, wave=wavelength)
osa.drawSpotDiagram()
sa = ShapeAnalysis(surf1.shape)
sa.plot(np.linspace(-1, 1, 10), np.linspace(-1, 1, 10), contours=100, ax=ax3)


plt.show()


plt.subplots_adjust(left=0.0, right=1.0, bottom=0.0, top=1.0)
# plots.drawLayout2d(ax2, s, [r2, r3, r4])
plt.subplots_adjust(left=0.0, right=1.0, bottom=0.0, top=1.0)
# plots.drawSpotDiagram(ax3, s, r3, -1)
# plt.subplots_adjust(left=0.0, right=1.0, bottom=0.0, top=1.0)
