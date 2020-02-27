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

import os
import math
import logging

from pprint import pprint

from distutils.version import StrictVersion

import numpy as np
from numpy import pi
import matplotlib.pyplot as plt
import matplotlib

from pyrateoptics.core.functionobject import FunctionObject
from pyrateoptics.core.base_ui import UIInterfaceClassWithOptimizableVariables
from pyrateoptics.core.serializer import Serializer, Deserializer
from pyrateoptics.optimize.optimize import Optimizer
from pyrateoptics.optimize.optimize_backends import ScipyBackend
from pyrateoptics.sampling2d.raster import RectGrid

from pyrateoptics import (build_rotationally_symmetric_optical_system,
                          listOptimizableVariables, draw)
from pyrateoptics.raytracer.globalconstants import canonical_ex, canonical_ey
from pyrateoptics.raytracer.globalconstants import Fline, dline, Cline
from pyrateoptics.raytracer.ray import RayBundle
from pyrateoptics.raytracer.analysis.ray_analysis import RayBundleAnalysis

logging.basicConfig(level=logging.INFO)

db_path = "refractiveindex.info-database/database"
serialization_path = "serialization_stuff/"

# drawing parameters
phi = 0.  # math.pi/4
pn = np.array([math.cos(phi), 0, math.sin(phi)])  # canonical_ex
up = canonical_ey

fig, axarr = plt.subplots(4, 1)
for ax in axarr:
    ax.axis('equal')
    if StrictVersion(matplotlib.__version__) < StrictVersion('2.0.0'):
        ax.set_axis_bgcolor('white')
    else:
        ax.set_facecolor('white')


# glass combinations from recipe:
# SK16 / BK7 / F5
# or
# LAK8 / BK7 / N-F2

# gcat = material_glasscat.refractiveindex_dot_info_glasscatalog(db_path)
# print gcat.findPagesWithLongNameContaining("BK7")


# Step 1: set up system of glass plates
########################################

rba = RayBundleAnalysis(None)

(s, seq) = build_rotationally_symmetric_optical_system(
        [(0,     0,     10.,    "N-SK16",         "lens1front", {}),
         (0,     0,     5,     None,             "lens1rear", {}),
         (0,     0,     5,     "N-BK7 (SCHOTT)", "elem2front", {}),
         (0,     0,     5,     "F5",             "elem2cement", {}),
         (0,     0,     5,     None,             "elem2rear", {}),
         (0,     0,     5,     None,             "stop", {"stop": True}),
         (0,     0,     5,     "F5",             "elem3front", {}),
         (0,     0,     5,     "N-BK7 (SCHOTT)", "elem3cement", {}),
         (0,     0,     5,     None,             "elem3rear", {}),
         (0,     0,     5,     "N-SK16",         "lens4front", {}),
         (0,     0,     5,     None,             "lens4rear", {}),
         (0,     0,     150.,     None,          "image", {})
         ], material_db_path=db_path, name="os")


def bundle_step1(nrays=100, rpup=7.5):
    """
    Creates an on-axis collimated RayBundle for step 1.
    """

    (px, py) = RectGrid().getGrid(nrays)
    o = np.vstack((rpup*px, rpup*py, np.zeros_like(px)))
    k = np.zeros_like(o)
    k[2, :] = 1.  # 2.*math.pi/wavelength
    E0 = np.cross(k, canonical_ex, axisa=0, axisb=0).T

    return RayBundle(o, k, E0, wave=dline)

initialbundle_step1 = bundle_step1()


# draw glass plates
s.draw2d(axarr[0], color="grey", vertices=50, plane_normal=pn, up=up)
r2 = s.seqtrace(bundle_step1(nrays=9), seq)
for r in r2:
    r.draw2d(axarr[0], color="blue", plane_normal=pn, up=up)


# Step 2: make last lens biconvex
#################################

# EFL = 1.5 * EFL of final objective = 150
# F/10 => EPD = 15
# monochrome d
# r1 = -r2

def meritfunction_step2(my_s):
    """
    Merit function (=function to be minimized) for step 2.
    """
    # initialbundle = bundle_step1()
    rpaths = my_s.seqtrace(initialbundle_step1, seq)
    x = rpaths[0].raybundles[-1].x[-1, 0, :]
    # y = rpaths[0].raybundles[-1].x[-1, 1, :]
    rba.raybundle = rpaths[0].raybundles[-1]

    merit = rba.getRMSspotSizeCentroid()  # np.sum(x**2 + y**2)
    merit += 1000000.*math.exp(-len(x))
    # TODO: adding the exp-x is a dirty trick.
    # The prefactor also has to be adapted for each system.
    # Is there a more elegant solution ?

    return merit


def updatefunction_allsteps(my_s):
    """
    Update coordinate systems of optical system.
    """
    my_s.rootcoordinatesystem.update()


# optimize
s.elements["stdelem"].surfaces["lens4front"].shape.curvature.toVariable()
s.elements["stdelem"].surfaces["lens4rear"].shape.curvature.toVariable()
s.elements["stdelem"].surfaces["lens4rear"].shape.curvature.\
    toPickup((FunctionObject("f = lambda x: -x"), "f"),
             (s.elements["stdelem"].surfaces["lens4front"].shape.curvature,))


listOptimizableVariables(s, max_line_width=80)

optimi = Optimizer(s, meritfunction_step2, backend=ScipyBackend(),
                   updatefunction=updatefunction_allsteps, name="step2")
s = optimi.run()

# draw system with last lens biconvex
s.draw2d(axarr[1], color="grey", vertices=50, plane_normal=pn, up=up)
r2 = s.seqtrace(bundle_step1(nrays=9), seq)  # trace again
for r in r2:
    r.draw2d(axarr[1], color="blue", plane_normal=pn, up=up)


# Step 3: colour, field, more variables
########################################


def bundles_step3(nrays=100, rpup=7.5, maxfield_deg=2.):
    """
    Creates an array of collimated RayBundle objects for step 3.
    """
    bundles = []
    fields_rad = np.array([0, 0.5*np.sqrt(2), 1]) * maxfield_deg * pi / 180.

    (px, py) = RectGrid().getGrid(nrays)

    for field in fields_rad:
        starty = -30 * np.tan(field)
        # TODO: this step explicitly sets the pupil
        # so the stop position is ignored. Better get Aimy to run ...

        o = np.vstack((rpup*px, rpup*py + starty, np.zeros_like(px)))
        k = np.zeros_like(o)
        k[1, :] = np.sin(field)
        k[2, :] = np.cos(field)
        E0 = np.cross(k, canonical_ex, axisa=0, axisb=0).T

        bundles += [RayBundle(o, k, E0, wave=Fline)]
        bundles += [RayBundle(o, k, E0, wave=dline)]
        bundles += [RayBundle(o, k, E0, wave=Cline)]

    return bundles


def final_meritfunction(my_s, rpup, fno, maxfield_deg):
    """
    Merit function (=function to be minimized) for step 3 and following.
    """
    merit = 0

    initialbundles = bundles_step3(rpup=rpup, maxfield_deg=maxfield_deg)
    for bundle in initialbundles:
        rpaths = my_s.seqtrace(bundle, seq)
        x = rpaths[0].raybundles[-1].x[-1, 0, :]
        rba.raybundle = rpaths[0].raybundles[-1]
        merit += rba.getRMSspotSizeCentroid()
        merit += 1./(len(x) + 1e-18)  # 10000.*math.exp(-len(x))

    # f number
    o = np.array([[0], [7.5], [0]])
    k = np.zeros_like(o)
    k[2, :] = 1
    E0 = np.cross(k, canonical_ex, axisa=0, axisb=0).T

    initial_marginal = RayBundle(o, k, E0, wave=dline)
    marginal_path = s.seqtrace(initial_marginal, seq)

    kyim = marginal_path[0].raybundles[-1].k[-1, 1, :]
    kzim = marginal_path[0].raybundles[-1].k[-1, 2, :]

    merit += 10 * sum(fno + .5 * np.real(kzim)/np.real(kyim))**2

    return merit


def meritfunction_step3(my_s):
    """
    Merit function for step 3
    """
    return final_meritfunction(my_s, rpup=7.5, fno=100/15., maxfield_deg=2.)


# optimize
s.elements["stdelem"].surfaces["lens1front"].shape.curvature.toVariable()
s.elements["stdelem"].surfaces["elem3front"].shape.curvature.toVariable()
s.elements["stdelem"].surfaces["elem3rear"].shape.curvature.toVariable()

# outer radii of both doulets should be symmetric
s.elements["stdelem"].surfaces["elem2front"].shape.curvature.\
    toPickup((FunctionObject("f = lambda x: -x"), "f"),
             (s.elements["stdelem"].surfaces["elem3rear"].shape.curvature,))

# inner radii of both doulets should be symmetric
s.elements["stdelem"].surfaces["elem2rear"].shape.curvature.\
    toPickup((FunctionObject("f = lambda x: -x"), "f"),
             (s.elements["stdelem"].surfaces["elem3front"].shape.curvature,))

s.elements["stdelem"].surfaces["image"].rootcoordinatesystem.decz.toVariable()


s.elements["stdelem"].surfaces["lens4front"].shape.curvature.toFixed()
s.elements["stdelem"].surfaces["lens4rear"].shape.curvature.toFixed()


optimi = Optimizer(s, meritfunction_step3, backend=ScipyBackend(),
                   updatefunction=updatefunction_allsteps, name="step3")
s = optimi.run()


s.draw2d(axarr[2], color="grey", vertices=50, plane_normal=pn, up=up)
b = bundles_step3(nrays=9)
for i in np.arange(len(b)):
    r2 = s.seqtrace(b[i], seq)
    color = ["blue", "green", "red"][int(i / 3)]
    for r in r2:
        r.draw2d(axarr[2], color=color, plane_normal=pn, up=up)


# Step 4: larger field
#######################

def meritfunction_step4a(my_s):
    """
    Merit function with increased field angle.
    """
    return final_meritfunction(my_s, rpup=7.5, fno=100/15., maxfield_deg=5.)


def meritfunction_step4b(my_s):
    """
    Merit function with most increased field angle.
    """
    return final_meritfunction(my_s, rpup=7.5, fno=100/15., maxfield_deg=7.)


s.elements["stdelem"].surfaces["lens4front"].shape.curvature.toVariable()
s.elements["stdelem"].surfaces["lens4rear"].shape.curvature.toVariable()

optimi = Optimizer(s, meritfunction_step4a, backend=ScipyBackend(),
                   updatefunction=updatefunction_allsteps, name="step4a")
s = optimi.run()
optimi = Optimizer(s, meritfunction_step4b, backend=ScipyBackend(),
                   updatefunction=updatefunction_allsteps, name="step4b")
s = optimi.run()


listOptimizableVariables(s, max_line_width=80)

# draw final result in overview plot fig and in fig2
fig2 = plt.figure()
ax2 = fig2.add_subplot(111)
if StrictVersion(matplotlib.__version__) < StrictVersion('2.0.0'):
    ax2.set_axis_bgcolor('white')
else:
    ax2.set_facecolor('white')
ax2.axis('equal')


s.draw2d(axarr[3], color="grey", vertices=50, plane_normal=pn, up=up)
s.draw2d(ax2,      color="grey", vertices=50, plane_normal=pn, up=up)
b = bundles_step3(nrays=9, rpup=7.5, maxfield_deg=5.)
for i in np.arange(len(b)):
    r2 = s.seqtrace(b[i], seq)
    color = ["blue", "green", "red"][int(i / 3)]
    for r in r2:
        r.draw2d(axarr[3], color=color, plane_normal=pn, up=up)

b = bundles_step3(nrays=36, rpup=7.5, maxfield_deg=5.)
r2 = s.seqtrace(b[0], seq)
for r in r2:
    r.draw2d(ax2, color="blue", plane_normal=pn, up=up)

plt.show()

system_gui_toplevel = UIInterfaceClassWithOptimizableVariables(
        s.elements["stdelem"].surfaces["elem2rear"].shape).query_for_dictionary()

pprint(system_gui_toplevel)


system_dump = Serializer(s)

try:
    os.mkdir(serialization_path)
except FileExistsError:
    pass

system_dump.save_json(serialization_path + "double_gauss.json")
system_dump.save_yaml(serialization_path + "double_gauss.yaml")

# TODO: later
#s2 = Deserializer.load_json(serialization_path + "double_gauss.json",
#                            True, True)
#draw(s2)
