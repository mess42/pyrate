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

######################################################
# Recipe "Linsensuppe a la Gauss" from Mo's Cookbook #
######################################################

# Step 0: Incredients
######################

import numpy as np
from numpy import pi
import matplotlib.pyplot as plt
import matplotlib
import math

from core.optical_system import OpticalSystem
from core.optical_element import OpticalElement
from core.surface import Surface
from core.surfShape import Conic
from core.localcoordinates import LocalCoordinates
from core.helpers import build_simple_optical_system
from core.helpers import collimated_bundle
from core.globalconstants import canonical_ex, canonical_ey
from core.ray import RayBundle
from core.optimize import Optimizer
from core.optimize_backends import ScipyBackend
from core.raster import RectGrid
from core.globalconstants import Fline, dline, Cline
from core import material_glasscat

from distutils.version import StrictVersion

db_path = "core/refractiveindex.info-database/database"

# drawing parameters
phi = 0.#math.pi/4
pn = np.array([math.cos(phi), 0, math.sin(phi)]) # canonical_ex
up = canonical_ey

fig, axarr = plt.subplots(4,1)
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

(s, seq) = build_simple_optical_system(
        [(0, 	0, 	10.,	"N-SK16", 		"lens1front"),
	 (0, 	0, 	5, 	None, 			"lens1rear"),
	 (0, 	0, 	5, 	"N-BK7 (SCHOTT)", 	"elem2front"),
	 (0, 	0, 	5, 	"F5", 			"elem2cement"),
	 (0, 	0, 	5, 	None, 			"elem2rear"),
	 (0, 	0, 	5, 	None, 			"stop"),
	 (0, 	0, 	5, 	"F5", 			"elem3front"),
	 (0, 	0, 	5, 	"N-BK7 (SCHOTT)", 	"elem3cement"),
	 (0, 	0, 	5, 	None, 			"elem3rear"),
         (0, 	0, 	5, 	"N-SK16", 		"lens4front"),
	 (0, 	0, 	5, 	None, 			"lens4rear"),
	 (0, 	0, 	150., 	None, 			"image")
         ], db_path)

def bundle_step1():
    """
    Creates an on-axis collimated RayBundle for step 1.
    """
    nrays = 100
    rpup = 7.5

    (px, py) = RectGrid().getGrid(nrays)
    o = np.vstack((rpup*px, rpup*py, np.zeros_like(px)))
    k = np.zeros_like(o)
    k[2,:] = 1. #2.*math.pi/wavelength
    E0 = np.cross(k, canonical_ex, axisa=0, axisb=0).T

    return RayBundle( o, k, E0, wave= dline)

# draw glass plates
s.draw2d(axarr[0], color="grey", vertices=50, plane_normal=pn, up=up)
r2 = s.seqtrace(bundle_step1(), seq)
for r in r2:
    r.draw2d(axarr[0], color="blue", plane_normal=pn, up=up) 



# Step 2: make last lens biconvex
##################################

# EFL = 1.5 * EFL of final objective = 150
# F/10 => EPD = 15
# monochrome d
# r1 = -r2

s.elements["stdelem"].surfaces["lens4front"].shape.curvature.changetype("variable")
s.elements["stdelem"].surfaces["lens4rear"].shape.curvature.changetype("variable")

def meritfunction_step2(s):
    """
    Merit function (=function to be minimized) for step 2.
    """
    initialbundle = bundle_step1()
    rpaths = s.seqtrace(initialbundle, seq)   
    x = rpaths[0].raybundles[-1].x[-1, 0, :]
    y = rpaths[0].raybundles[-1].x[-1, 1, :]
    
    merit  = np.sum(x**2 + y**2) + 1000000.*math.exp(-len(x))
    # TODO: adding the exp-x is a dirty trick. The prefactor also has to be adapted for each system. Is there a more elegant solution ?

    # biconvex lens with same radii
    merit += 10000* (  s.elements["stdelem"].surfaces["lens4front"].shape.curvature.parameters["value"] \
                     + s.elements["stdelem"].surfaces["lens4rear"].shape.curvature.parameters["value"]  )**2
    # TODO: OptimizableVariable() documentation does not explain what an OptimizableVariable is and makes it cumbersome to use

    return merit

# optimize
optimi = Optimizer(s, meritfunction_step2, backend=ScipyBackend(), updatefunction=s.rootcoordinatesystem.update()) 
# TODO: Optimizer() is not well documented
s = optimi.run()

# draw system with last lens biconvex
s.draw2d(axarr[1], color="grey", vertices=50, plane_normal=pn, up=up) # try for phi=0.
r2 = s.seqtrace(bundle_step1(), seq) # trace again
for r in r2:
    r.draw2d(axarr[1], color="blue", plane_normal=pn, up=up) 



# Step 3: colour, field, more variables
########################################

def bundles_step3(nrays = 100, rpup = 7.5, maxfield_deg=2.):
    """
    Creates an array of collimated RayBundle objects for step 3.
    """
    bundles = []
    fields_rad = np.array([0,.5*np.sqrt(2),1]) * maxfield_deg * pi /180.

    (px, py) = RectGrid().getGrid(nrays)

    for field in fields_rad:
        starty = 20 * np.tan(field) # TODO: this step explicitly sets the pupil, so the stop position is ignored. Better get Aimy to run ...
        o = np.vstack((rpup*px, rpup*py + starty, np.zeros_like(px)))
        k = np.zeros_like(o)
        k[1,:] = np.sin( field )
        k[2,:] = np.cos( field )
        E0 = np.cross(k, canonical_ex, axisa=0, axisb=0).T 
        
        bundles += [RayBundle( o, k, E0, wave= Fline)]
        bundles += [RayBundle( o, k, E0, wave= dline)]
        bundles += [RayBundle( o, k, E0, wave= Cline)]

    return bundles




s.draw2d(axarr[2], color="grey", vertices=50, plane_normal=pn, up=up) # try for phi=0.
for b in bundles_step3():
    r2 = s.seqtrace(b, seq) # trace again
    for r in r2:
        r.draw2d(axarr[2], color="blue", plane_normal=pn, up=up) 


plt.show()



