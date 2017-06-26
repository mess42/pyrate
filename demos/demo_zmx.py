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

import numpy as np
import math
import sys

from core.localcoordinates import LocalCoordinates
from core.material_isotropic import ConstantIndexGlass
from core.globalconstants import standard_wavelength, numerical_tolerance
from core.ray import RayBundle
from core.zmx import ZMXParser

from core import raster

import matplotlib
import matplotlib.pyplot as plt
from distutils.version import StrictVersion


# download ZMX files from e.g.:
# http://astro.dur.ac.uk/~rsharp/opticaldesign/
# some good demonstration of coordinate breaks is: FIELDROTATOR-LECT5.ZMX
print(len(sys.argv))
if len(sys.argv) != 5:
    print("usage:")
    print("python -m demos.demo_zmx file.zmx u|a entrance_pupil_diameter num_rays_for_yfan")
    print("only collimated light")
    print("file.zmx relative to local directory")
    print("u for utf16 file; a for ascii file (try both)")
    exit()

file_to_read = sys.argv[1]
u_or_a = sys.argv[2]
enpd = float(sys.argv[3])
num_rays = int(sys.argv[4])

ascii_read = True if u_or_a == "a" else False

p = ZMXParser(file_to_read, ascii=ascii_read)
lctmp = LocalCoordinates("tmp")

#matdict = {}
matdict = {"BK7":ConstantIndexGlass(lctmp, 1.5168)}
#matdict = {"LAFN21":ConstantIndexGlass(lctmp, 1.788), "SF53":ConstantIndexGlass(lctmp, 1.72)}    

(s, seq) = p.createOpticalSystem(matdict)


rstobj = raster.MeridionalFan()
(px, py) = rstobj.getGrid(num_rays)

rpup = enpd*0.5 #7.5
o = np.vstack((rpup*px, rpup*py, -5.*np.ones_like(px)))

k = np.zeros_like(o)
k[1,:] = 2.*math.pi/standard_wavelength*math.sin(0.0)
k[2,:] = 2.*math.pi/standard_wavelength*math.cos(0.0)

ey = np.zeros_like(o)
ey[1,:] =  1.

E0 = np.cross(k, ey, axisa=0, axisb=0).T

initialbundle = RayBundle(x0=o, k0=k, Efield0=E0, wave=standard_wavelength)
rays = s.seqtrace(initialbundle, seq)




fig = plt.figure(1)
ax = fig.add_subplot(111)

ax.axis('equal')
if StrictVersion(matplotlib.__version__) < StrictVersion('2.0.0'):
    ax.set_axis_bgcolor('white')
else:
    ax.set_facecolor('white')

for r in rays:
    r.draw2d(ax, color="blue")

s.draw2d(ax, color="grey", vertices=50, inyzplane=False)

plt.show()