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
import matplotlib.pyplot as plt

import matplotlib
from distutils.version import StrictVersion


from pyrateoptics.core import raster
from pyrateoptics.core.material_isotropic import ConstantIndexGlass
from pyrateoptics.core.surfShape import Conic, Asphere
from pyrateoptics.core.optical_element import OpticalElement
from pyrateoptics.core.surface import Surface
from pyrateoptics.core.optical_system import OpticalSystem
from pyrateoptics.core.ray import RayBundle

from pyrateoptics.core.aperture import CircularAperture
from pyrateoptics.core.localcoordinates import LocalCoordinates

from pyrateoptics.core.globalconstants import canonical_ey

from pyrateoptics.core.optimize import Optimizer
from pyrateoptics.core.optimize_backends import ScipyBackend

import math
import logging
logging.basicConfig(level=logging.DEBUG)

wavelength = 0.5876e-3

# definition of optical system
s = OpticalSystem() 

lc0 = s.addLocalCoordinateSystem(LocalCoordinates(name="stop", decz=0.0), refname=s.rootcoordinatesystem.name)
lc1 = s.addLocalCoordinateSystem(LocalCoordinates(name="surf1", decz=5.0), refname=lc0.name) # objectDist
lc2 = s.addLocalCoordinateSystem(LocalCoordinates(name="surf2", decz=20.0), refname=lc1.name)
lc3 = s.addLocalCoordinateSystem(LocalCoordinates(name="image", decz=100.0), refname=lc2.name)


stopsurf = Surface(lc0)
frontsurf = Surface(lc1, shape=Conic(lc1), apert=CircularAperture(lc1, 12.7))
backsurf = Surface(lc2, shape=Asphere(lc2, curv=-1./50.0, cc=-1., coefficients=[0.0, 0.0, 0.0]), apert=CircularAperture(lc2, 12.7))
image = Surface(lc3)


elem = OpticalElement(lc0, name="asphereelement")

bk7 = ConstantIndexGlass(lc1, n=1.5168)

elem.addMaterial("BK7", bk7)

elem.addSurface("stop", stopsurf, (None, None))
elem.addSurface("front", frontsurf, (None, "BK7"))
elem.addSurface("rear", backsurf, ("BK7", None))
elem.addSurface("image", image, (None, None))

s.addElement("asph", elem)

rstobj = raster.MeridionalFan()
(px, py) = rstobj.getGrid(20)

rpup = 11.43
o = np.vstack((rpup*px, rpup*py, -5.*np.ones_like(px)))

k = np.zeros_like(o)
k[2,:] = 1. #2.*math.pi/wavelength

ey = np.zeros_like(o)
ey[1,:] =  1.

E0 = np.cross(k, ey, axisa=0, axisb=0).T

sysseq = [("asph", [("stop", {"is_stop":True}), 
                    ("front", {}), 
                    ("rear", {}), 
                    ("image", {})])]

phi = 5.*math.pi/180.0

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

backsurf.shape.params["curv"].changetype("variable")
backsurf.shape.params["cc"].changetype("variable")
# A2 not variable
backsurf.shape.params["A4"].changetype("variable")
backsurf.shape.params["A6"].changetype("variable")

opt_backend = ScipyBackend(method='Nelder-Mead', tol=1e-9)
optimi = Optimizer(s, meritfunctionrms, opt_backend, name="Nelder-Mead Optimizer")
s = optimi.run()

r2 = s.seqtrace(initialbundle, sysseq)


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

for r in r2:
    r.draw2d(ax, color="blue", plane_normal=pn, up=up) 
for e in s.elements.itervalues():
    for surfs in e.surfaces.itervalues():
        surfs.draw2d(ax, color="grey", vertices=50, plane_normal=pn, up=up) # try for phi=0.
        #surfs.draw2d(ax, color="grey", inyzplane=False, vertices=50, plane_normal=pn, up=up) # try for phi=pi/4


plt.show()


