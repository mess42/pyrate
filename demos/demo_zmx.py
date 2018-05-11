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
logging.basicConfig(level=logging.DEBUG)

from pyrateoptics.raytracer.localcoordinates import LocalCoordinates
from pyrateoptics.material.material_isotropic import ConstantIndexGlass
from pyrateoptics.raytracer.globalconstants import standard_wavelength
from pyrateoptics.raytracer.ray import RayBundle
from pyrateoptics.io.zmx import ZMXParser

from pyrateoptics.sampling2d import raster
from pyrateoptics import collimated_bundle, draw



# download ZMX files from e.g.:
# http://astro.dur.ac.uk/~rsharp/opticaldesign/
# some good demonstration of coordinate breaks is: FIELDROTATOR-LECT5.ZMX

if len(sys.argv) != 4:
    print("usage:")
    print("python -m demos.demo_zmx file.zmx entrance_pupil_diameter num_rays_for_yfan")
    print("only collimated light")
    print("file.zmx relative to local directory")
    exit()

file_to_read = sys.argv[1]
enpd = float(sys.argv[2])
num_rays = int(sys.argv[3])

p = ZMXParser(file_to_read, name='ZMXParser')
lctmp = LocalCoordinates("tmp")

matdict = {}
#matdict = {"BK7":ConstantIndexGlass(lctmp, 1.5168)}
#matdict = {"LAFN21":ConstantIndexGlass(lctmp, 1.788), "SF53":ConstantIndexGlass(lctmp, 1.72)}    

(s, seq) = p.createOpticalSystem(matdict)
#(o, k, E0) = collimated_bundle(11, {"opticalsystem":s, "radius":enpd*0.5, "startz":-5., "raster":raster.MeridionalFan()}, wave=standard_wavelength)

initialbundles_dict = p.createInitialBundle()
print(seq)
ray_paths = []
for d in initialbundles_dict:
    d["opticalsystem"] = s
    d["raster"] = raster.MeridionalFan()
    (o, k, E0) = collimated_bundle(11, d, wave=standard_wavelength)
    initialbundle = RayBundle(x0=o, k0=k, Efield0=E0, wave=standard_wavelength)
    ray_paths.append(s.seqtrace(initialbundle, seq))

draw(s, ray_paths)
