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

# Incredients:
import numpy as np
from core.optical_system import OpticalSystem
from core.optical_element import OpticalElement
from core.surface import Surface
from core.surfShape import Conic
from core.localcoordinates import LocalCoordinates
from core.helpers import build_simple_optical_system
from core import material_glasscat

db_path = "core/refractiveindex.info-database/database"

# pre-design: find suitable glasses
# gcat = material_glasscat.refractiveindex_dot_info_glasscatalog(db_path) 
# print gcat.findPagesWithLongNameContaining("BK7")


# step 0: set up system of glass plates
(s, seq) = build_simple_optical_system(
        [(0, 	0, 	1000., 	"N-SK16", 		"lens1front"),
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
	 (0, 	0, 	100., 	None, 			"image")
         ], db_path)

