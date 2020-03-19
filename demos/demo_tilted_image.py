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

import numpy as np

from pyrateoptics import build_rotationally_symmetric_optical_system, draw
from pyrateoptics.raytracer.globalconstants import degree
from pyrateoptics.raytracer.analysis.optical_system_analysis import\
    OpticalSystemAnalysis
from pyrateoptics.raytracer.ray import RayBundle
from pyrateoptics.sampling2d.raster import MeridionalFan
from pyrateoptics.raytracer.helpers import build_pilotbundle_complex

alpha = 10.*degree

epd = 5.

(s, seq) = build_rotationally_symmetric_optical_system(
        [(0, 0, 0., None, "object", {}),
         (100., 0, 5, 1.5, "lens1front", {"is_stop": True}),
         (0., 0, 5, None, "lens1rear", {}),
         (0, 0, 196.228, None, "image", {})], name="os")

imsurf = s.elements["stdelem"].surfaces["image"]
imsurf.rootcoordinatesystem.tiltx.set_value(alpha)
imsurf.rootcoordinatesystem.update()
objsurf = s.elements["stdelem"].surfaces["object"]
osa = OpticalSystemAnalysis(s, seq, name="Analysis")
(x01, k01, E01) = osa.collimated_bundle(11, {"radius": epd,
                                             "raster": MeridionalFan()})
(x02, k02, E02) = osa.collimated_bundle(11, {"radius": epd,
                                             "raster": MeridionalFan(),
                                             "anglex": 1.*degree})
(x03, k03, E03) = osa.collimated_bundle(11, {"radius": epd,
                                             "raster": MeridionalFan(),
                                             "anglex": -1.*degree})

mybundle1 = RayBundle(x01, k01, E01)
mybundle2 = RayBundle(x02, k02, E02)
mybundle3 = RayBundle(x03, k03, E03)

raypaths1 = s.seqtrace(mybundle1, seq)
raypaths2 = s.seqtrace(mybundle2, seq)
raypaths3 = s.seqtrace(mybundle3, seq)

obj_dx = 0.1
obj_dphi = 0.1*degree

pilotbundles = build_pilotbundle_complex(objsurf, s.material_background,
                                         (obj_dx, obj_dx),
                                         (obj_dphi, obj_dphi),
                                         num_sampling_points=3)
(m_obj_stop, m_stop_img) = s.extractXYUV(pilotbundles[-1], seq)

print(np.array_str(m_obj_stop, precision=5, suppress_small=True))
print(np.array_str(m_stop_img, precision=5, suppress_small=True))

draw(s, [raypaths1, raypaths2, raypaths3])
