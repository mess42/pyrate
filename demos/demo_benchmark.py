#!/usr/bin/env/python
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

import time
import sys
import logging

from pyrateoptics import build_rotationally_symmetric_optical_system, draw
from pyrateoptics.sampling2d import raster
from pyrateoptics.raytracer.ray import RayBundle
from pyrateoptics.raytracer.globalconstants import standard_wavelength, degree
from pyrateoptics.raytracer.analysis.optical_system_analysis import\
    OpticalSystemAnalysis

from pyrateoptics.raytracer.config import ConfigFile

logging.basicConfig(level=logging.INFO)

wavelength = standard_wavelength

# definition of optical system

(s, seq) = build_rotationally_symmetric_optical_system(
        [(-5.922, 0, 2.0, 1.7, "surf1", {}),
         (-3.160, 0, 3.0, None, "surf2", {}),
         (15.884, 0, 5.0, 1.7, "surf3", {}),
         (-12.756, 0, 3.0, None, "surf4", {}),
         (0, 0, 3.0, None, "stop", {"is_stop": True}),
         (3.125, 0, 2.0, 1.5, "surf5", {}),
         (1.479, 0, 3.0, None, "surf6", {}),
         (0, 0, 19.0, None, "surf7", {})
         ], material_db_path=ConfigFile().get_refractive_index_database_path())


def mytiming():
    if sys.version_info.major >= 3:
        return time.perf_counter()
    else:
        return time.clock()


nrays = 100000
nrays_draw = 21

osa = OpticalSystemAnalysis(s, seq, name="Analysis")

(x0, k0, E0) = osa.divergent_bundle(nrays,
                                    {"radius": 10.*degree,
                                     "raster": raster.RectGrid()})
t0 = mytiming()
initialraybundle = RayBundle(x0=x0, k0=k0, Efield0=E0)
t1 = mytiming()
raypath = s.seqtrace(initialraybundle, seq)
t2 = mytiming()
logging.info("benchmark : " + str(t2 - t1) +
             " s for tracing " + str(nrays) + " rays through " +
             str(len(s.elements["stdelem"].surfaces) - 1) + " surfaces.")
logging.info("That is " +
             str(int(round(nrays * (len(s.elements["stdelem"].surfaces) - 1)
                           / (t2 - t1)))) +
             " ray-surface-operations per second")

# plot

(x0_draw, k0_draw, E0_draw) =\
    osa.divergent_bundle(nrays_draw,
                         {"radius": 10.*degree,
                          "raster": raster.MeridionalFan()})
initialraybundle_draw = RayBundle(x0=x0_draw, k0=k0_draw, Efield0=E0_draw)
raypath_draw = s.seqtrace(initialraybundle_draw, seq)


draw(s, raypath_draw)
