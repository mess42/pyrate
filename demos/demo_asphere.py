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

import numpy as np


from pyrateoptics.raytracer.ray import RayBundle
from pyrateoptics.analysis.optical_system_analysis import OpticalSystemAnalysis

from pyrateoptics.optimize.optimize import Optimizer
from pyrateoptics.optimize.optimize_backends import ScipyBackend

from pyrateoptics import build_simple_optical_system, draw, raytrace

import logging
logging.basicConfig(level=logging.DEBUG)

wavelength = 0.5876e-3

# definition of optical system
(s, sysseq) = build_simple_optical_system(
                [
                    ({"shape": "Conic"}, {"decz":0.0}, None, "stop", {"is_stop":True}),
                    ({"shape": "Conic"}, {"decz":5.0}, 1.5168, "front", {}),
                    ({"shape": "Asphere", "curv": -1./50., 
                          "cc": -1., "coefficients": [0.0, 0.0, 0.0]},
                            {"decz":20.0}, None, "back", {}),
                    ({"shape": "Conic"}, {"decz":100.0}, None, "image", {})
                ],
                )

osa = OpticalSystemAnalysis(s, sysseq, name="Analysis")

(o, k, E0) = osa.collimated_bundle(121, {"startz":-5., "radius":11.43}, wave=wavelength)
initialbundle = RayBundle(x0=o, k0=k, Efield0=E0, wave=wavelength)

#initialbundle = generatebundle(openangle=10.*math.pi/180, numrays=121)

def meritfunctionrms(s):
    initialbundle_local = RayBundle(x0=o, k0=k, Efield0=E0, wave=wavelength)
    rpaths = s.seqtrace(initialbundle_local, sysseq)
    # other constructions lead to fill up of initial bundle with intersection values
    
    # for glassy asphere only one path necessary
    x = rpaths[0].raybundles[-1].x[-1, 0, :]
    y = rpaths[0].raybundles[-1].x[-1, 1, :]
    
    res = np.sum(x**2 + y**2)
    
    return res

backsurf = s.elements["stdelem"].surfaces["back"]
backsurf.shape.params["curv"].changetype("variable")
backsurf.shape.params["cc"].changetype("variable")
# A2 not variable
backsurf.shape.params["A4"].changetype("variable")
backsurf.shape.params["A6"].changetype("variable")

opt_backend = ScipyBackend(method='Nelder-Mead', tol=1e-9)
optimi = Optimizer(s, meritfunctionrms, opt_backend, name="Nelder-Mead Optimizer")
s = optimi.run()

r2 = s.seqtrace(initialbundle, sysseq)

draw(s, r2)


