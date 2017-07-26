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
import matplotlib.pyplot as plt
import matplotlib
import math

from core.optical_system import OpticalSystem
from core.optical_element import OpticalElement
from core.surface import Surface
from core.surfShape import Conic
from core.localcoordinates import LocalCoordinates
from core.helpers import build_simple_optical_system
from core.globalconstants import canonical_ey
from core import material_glasscat

from distutils.version import StrictVersion

db_path = "core/refractiveindex.info-database/database"

# glass combinations from recipe:
# SK16 / BK7 / F5
# or
# LAK8 / BK7 / N-F2

# gcat = material_glasscat.refractiveindex_dot_info_glasscatalog(db_path) 
# print gcat.findPagesWithLongNameContaining("BK7")




# Step 1: set up system of glass plates
########################################

(s, seq) = build_simple_optical_system(
        [(0, 	0, 	10000.,	"N-SK16", 		"lens1front"),
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


# helper.collimated_bundle









fig = plt.figure(1)
ax = fig.add_subplot(111)
ax.axis('equal')
if StrictVersion(matplotlib.__version__) < StrictVersion('2.0.0'):
    ax.set_axis_bgcolor('white')
else:
    ax.set_facecolor('white')


phi = 0.#math.pi/4
pn = np.array([math.cos(phi), 0, math.sin(phi)]) # canonical_ex
up = canonical_ey

#for r in raypath_draw:
#    r.draw2d(ax, color="blue", plane_normal=pn, up=up) 

s.draw2d(ax, color="grey", vertices=50, plane_normal=pn, up=up) # try for phi=0.

plt.show()



