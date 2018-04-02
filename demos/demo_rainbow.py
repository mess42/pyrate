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
from pyrateoptics.core.material_glasscat import refractiveindex_dot_info_glasscatalog, CatalogMaterial
from pyrateoptics.core.surfShape import Asphere
from pyrateoptics.core.optical_element import OpticalElement
from pyrateoptics.core.surface import Surface
from pyrateoptics.core.optical_system import OpticalSystem
from pyrateoptics.core.ray import RayBundle

from pyrateoptics.core.aperture import CircularAperture
from pyrateoptics.core.localcoordinates import LocalCoordinates

from pyrateoptics.core.globalconstants import canonical_ey, degree

import math
import logging
logging.basicConfig(level=logging.DEBUG)

wavelength = 0.5876e-3

wave_red = 0.700e-3
wave_blue = 0.470e-3




# definition of optical system
s = OpticalSystem() 

dropletradius = 0.1

lc0 = s.addLocalCoordinateSystem(LocalCoordinates(name="stop", decz=0.0), refname=s.rootcoordinatesystem.name)
lccomprism = s.addLocalCoordinateSystem(LocalCoordinates(name="dropletcenter", decz=2.*dropletradius), refname=lc0.name)

lc1 = s.addLocalCoordinateSystem(LocalCoordinates(name="surf1", decz=-dropletradius), refname=lccomprism.name) # objectDist
lc2 = s.addLocalCoordinateSystem(LocalCoordinates(name="surf2", decz=dropletradius), refname=lccomprism.name)
lc3 = s.addLocalCoordinateSystem(LocalCoordinates(name="surf3", decz=0), refname=lccomprism.name)
lc4 = s.addLocalCoordinateSystem(LocalCoordinates(name="image", decz=-2.*dropletradius), refname=lccomprism.name)


stopsurf = Surface(lc0, apert=CircularAperture(lc0, 7*dropletradius))
frontsurf = Surface(lc1, shape=Asphere(lc1, curv=1./dropletradius), apert=CircularAperture(lc1, dropletradius))
rearsurf = Surface(lc2, shape=Asphere(lc2, curv=-1./dropletradius), apert=CircularAperture(lc2, dropletradius))
midsurf = Surface(lc3, shape=Asphere(lc3, curv=0), apert=CircularAperture(lc3, dropletradius))

image = Surface(lc4, apert=CircularAperture(lc4, 7.*dropletradius))


elem = OpticalElement(lc0, name="droplet")



try:
    database_basepath = "refractiveindex.info-database/database"     
    shelf = "3d"
    book  = "liquids"
    page  = "water"    


    gcat = refractiveindex_dot_info_glasscatalog(database_basepath)

    waterdict = gcat.getMaterialDict(shelf, book, page)
    water = CatalogMaterial(lc0, waterdict, name="water (catalogue)")

except:
    logging.warn("refractive index database not found. please download it and symlink\n\
            to it in your local pyrate directory")
    water = ConstantIndexGlass(lc0, n=1.336, name="water (failsafe)")

logging.info("wavelength %f, index %f" % (wave_red, water.getIndex(None, wave_red).real))
logging.info("wavelength %f, index %f" % (wave_blue, water.getIndex(None, wave_blue).real))

elem.addMaterial("water", water)

elem.addSurface("stop", stopsurf, (None, None))
elem.addSurface("surf1", frontsurf, (None, "water"))
elem.addSurface("surf2", rearsurf, ("water", "water"))
elem.addSurface("surf4", frontsurf, ("water", None))
elem.addSurface("image", image, (None, None))

s.addElement("droplet", elem)

rstobj = raster.MeridionalFan()
(px, py) = rstobj.getGrid(11)

rpup = dropletradius*0.05
oy = dropletradius*0.9
o = np.vstack((rpup*px, rpup*py + oy, 0*np.ones_like(px)))

kangle = -12.*degree

kwave_red = 1.0
k_red = np.zeros_like(o)
k_red[1,:] = kwave_red*math.sin(kangle)
k_red[2,:] = kwave_red*math.cos(kangle)

kwave_blue = 1.0
k_blue = np.zeros_like(o)
k_blue[1,:] = kwave_blue*math.sin(kangle)
k_blue[2,:] = kwave_blue*math.cos(kangle)


ey = np.zeros_like(o)
ey[1,:] =  1.

E0_red = np.cross(k_red, ey, axisa=0, axisb=0).T
E0_blue = np.cross(k_blue, ey, axisa=0, axisb=0).T

sysseq = [("droplet", 
           [
                ("stop", {"is_stop":True}), 
                ("surf1", {}), 
#                ("surf3", {}), 
                ("surf2", {"is_mirror":True}),
#                ("surf3", {}), 
                ("surf4", {}), 
                ("image", {})])]

sysseq2nd = [("droplet", [
                ("stop", {"is_stop":True}), 
                ("surf1", {}), 
#                ("surf3", {}), 
                ("surf2", {"is_mirror":True}), 
#                ("surf3", {}), 
                ("surf4", {}), 
                ("image", {})])]


initialbundle_red = RayBundle(x0=o, k0=k_red, Efield0=E0_red, wave=wave_red)
initialbundle_blue = RayBundle(x0=o, k0=k_blue, Efield0=E0_blue, wave=wave_blue)
r_red = s.seqtrace(initialbundle_red, sysseq)
r_blue = s.seqtrace(initialbundle_blue, sysseq)

#r_red2nd = s.seqtrace(initialbundle_red, sysseq2nd)
#r_blue2nd = s.seqtrace(initialbundle_blue, sysseq2nd)


fig = plt.figure(1)
ax = fig.add_subplot(111)

ax.axis('equal')
if StrictVersion(matplotlib.__version__) < StrictVersion('2.0.0'):
    ax.set_axis_bgcolor('white')
else:
    ax.set_facecolor('white')


phi = 0 #math.pi/4
pn = np.array([math.cos(phi), 0, math.sin(phi)]) # canonical_ex
up = canonical_ey

for r in r_red: 
    r.draw2d(ax, color="red", plane_normal=pn, up=up)
for r in r_blue:
    r.draw2d(ax, color="blue", plane_normal=pn, up=up) 

#for r in r_red2nd:
#    r.draw2d(ax, color="red", plane_normal=pn, up=up) 
#for r in r_blue2nd:
#    r.draw2d(ax, color="blue", plane_normal=pn, up=up) 


s.draw2d(ax, color="grey", vertices=50, plane_normal=pn, up=up) # try for phi=0.
#s.draw2d(ax, color="grey", inyzplane=False, vertices=50, plane_normal=pn, up=up) # try for phi=pi/4


plt.show()


