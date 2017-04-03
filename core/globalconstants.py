#!/usr/bin/env/python
"""
Pyrate - Optical raytracing based on Python

Copyright (C) 2014 Moritz Esslinger moritz.esslinger@web.de
               and Johannes Hartung j.hartung@gmx.net
               and    Uwe Lippmann  uwe.lippmann@web.de

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
import numpy as np

standard_wavelength = 0.5876e-3 # (d line) in mm

c0 = 299792458 # definition, m/s
mu0 = 4.*math.pi*1e-7 # N/A**2
eps0 = 1./(c0**2*mu0) # Vs/(Am)

canonical_ex = np.array([1, 0, 0])
canonical_ey = np.array([0, 1, 0])
canonical_ez = np.array([0, 0, 1])
