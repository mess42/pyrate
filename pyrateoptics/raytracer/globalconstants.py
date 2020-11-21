#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Pyrate - Optical raytracing based on Python

Copyright (C) 2014-2020
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

import math
import numpy as np

# spectral lines
iline = 0.3650E-3  # i  (mm)
hline = 0.4047E-3  # h  (mm)
gline = 0.4358E-3  # g  (mm)
Fprimeline = 0.4800E-3  # F' (mm)
Fline = 0.4861E-3  # F  (mm)
eline = 0.5461E-3  # e  (mm)
dline = 0.5876E-3  # d  (mm)
Dline = 0.5893E-3  # D  (mm)
Cprimeline = 0.6438E-3  # C' (mm)
Cline = 0.6563E-3  # C  (mm)
rline = 0.7065E-3  # r  (mm)
sline = 0.8521E-3  # s  (mm)
tline = 1.0140E-3  # t  (mm)

standard_wavelength = dline

c0 = 299792458  # definition, m/s
mu0 = 4.*math.pi*1e-7  # N/A**2
eps0 = 1./(c0**2*mu0)  # Vs/(Am)

canonical_ex = np.array([1, 0, 0])
canonical_ey = np.array([0, 1, 0])
canonical_ez = np.array([0, 0, 1])

degree = math.pi/180.0

numerical_tolerance = 1e-17
