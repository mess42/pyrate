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
from pyrateoptics import draw

import matplotlib.pyplot as plt

from pyrateoptics.analysis.optical_system_analysis import OpticalSystemAnalysis

import argparse

# download ZMX files from e.g.:
# http://astro.dur.ac.uk/~rsharp/opticaldesign/
# some good demonstration of coordinate breaks is: FIELDROTATOR-LECT5.ZMX

parser = argparse.ArgumentParser(description="Read ZMX files.")
parser.add_argument("file", help="File to be interpreted", type=str)
parser.add_argument("--bundletype", nargs='?', help="Bundle type", type=str, default='collimated')
parser.add_argument("--epd", nargs='?', help="Entrance pupil diameter", type=float, default=1.0)
parser.add_argument("--numrays", nargs='?', help="Number of rays", type=int, default=11)
parser.add_argument("--showspot", help="Show spot diagram?", action="store_true")
parser.add_argument("--anglex", help="Angle", type=float, default=0.0)
parsed = parser.parse_args()

# TODO: add materials via command line

show_spot = parsed.showspot
file_to_read = parsed.file
enpd = parsed.epd
num_rays = parsed.numrays
bundletype = parsed.bundletype
anglex = parsed.anglex

p = ZMXParser(file_to_read, name='ZMXParser')
lctmp = LocalCoordinates("tmp")

matdict = {"BK7":ConstantIndexGlass(lctmp, 1.5168), 
           "LAFN21":ConstantIndexGlass(lctmp, 1.788),
           "SF53":ConstantIndexGlass(lctmp, 1.72)}    

(s, seq) = p.createOpticalSystem(matdict)

if s is None:
    sys.exit()

initialbundles_dict = [p.createInitialBundle("N")[0]]

osa = OpticalSystemAnalysis(s, seq, name="Analysis")

ray_paths = []

if initialbundles_dict == []:
    initialbundles_dict = [{"radius":enpd*0.5}]

for d in initialbundles_dict:
    if show_spot:
        d["raster"] = raster.RectGrid()
        osa.aim(num_rays*num_rays, d, wave=standard_wavelength)
        osa.drawSpotDiagram()
    else:
        d["raster"] = raster.MeridionalFan()
        d["anglex"] = anglex        
        osa.aim(num_rays, d, wave=standard_wavelength)
        ray_paths.append(osa.trace()[0])

if not show_spot:
    draw(s, ray_paths)
else:
    plt.show()
osa.prettyprint()

