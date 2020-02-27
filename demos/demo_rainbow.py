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
import logging
import json
import yaml

from pprint import pprint


from pyrateoptics.sampling2d import raster
from pyrateoptics.raytracer.material.material_isotropic import\
    ConstantIndexGlass
from pyrateoptics.raytracer.material.material_glasscat import\
    refractiveindex_dot_info_glasscatalog, CatalogMaterial
from pyrateoptics.raytracer.surface_shape import Asphere
from pyrateoptics.raytracer.optical_element import OpticalElement
from pyrateoptics.raytracer.surface import Surface
from pyrateoptics.raytracer.optical_system import OpticalSystem

from pyrateoptics.raytracer.aperture import CircularAperture
from pyrateoptics.raytracer.localcoordinates import LocalCoordinates

from pyrateoptics.core.base_ui import UIInterfaceClassWithOptimizableVariables
from pyrateoptics.core.serializer import Serializer


from pyrateoptics.raytracer.globalconstants import degree
from pyrateoptics import raytrace, draw

logging.basicConfig(level=logging.DEBUG)

wavelength = 0.5876e-3

wave_red = 0.700e-3
wave_blue = 0.470e-3

# definition of optical system
s = OpticalSystem.p()

dropletradius = 0.1

lc0 = s.addLocalCoordinateSystem(
            LocalCoordinates.p(name="stop", decz=0.0),
            refname=s.rootcoordinatesystem.name)
lccomprism = s.addLocalCoordinateSystem(
            LocalCoordinates.p(name="dropletcenter", decz=2.*dropletradius),
            refname=lc0.name)

lc1 = s.addLocalCoordinateSystem(LocalCoordinates.p(name="surf1",
                                                    decz=-dropletradius),
                                 refname=lccomprism.name)  # objectDist
lc2 = s.addLocalCoordinateSystem(
            LocalCoordinates.p(name="surf2", decz=dropletradius),
            refname=lccomprism.name)
lc3 = s.addLocalCoordinateSystem(
            LocalCoordinates.p(name="surf3", decz=0),
            refname=lccomprism.name)
lc4 = s.addLocalCoordinateSystem(
            LocalCoordinates.p(name="image", decz=-2.*dropletradius),
            refname=lccomprism.name)


stopsurf = Surface.p(lc0,
                     aperture=CircularAperture.p(lc0, maxradius=7*dropletradius))
frontsurf = Surface.p(lc1, shape=Asphere.p(lc1, curv=1./dropletradius),
                      aperture=CircularAperture.p(lc1, maxradius=dropletradius))
rearsurf = Surface.p(lc2, shape=Asphere.p(lc2, curv=-1./dropletradius),
                     aperture=CircularAperture.p(lc2, maxradius=dropletradius))
midsurf = Surface.p(lc3, shape=Asphere.p(lc3, curv=0),
                    aperture=CircularAperture.p(lc3, maxradius=dropletradius))

image = Surface.p(lc4,
                aperture=CircularAperture.p(lc4, maxradius=7.*dropletradius))


elem = OpticalElement.p(lc0, name="droplet")


database_basepath = "refractiveindex.info-database/database"
shelf = "3d"
book = "liquids"
page = "water"

try:
    gcat = refractiveindex_dot_info_glasscatalog(database_basepath)
    waterdict = gcat.getMaterialDict(shelf, book, page)
    water = CatalogMaterial.p(lc0, waterdict, name="water (catalogue)")
except KeyError:
    logging.warning("refractive index database not found. please download it\
                     and symlink\nto it in your local pyrate directory")
    water = ConstantIndexGlass.p(lc0, n=1.336, name="water (failsafe)")

logging.info("wavelength %f, index %f" %
             (wave_red, water.getIndex(None, wave_red).real))
logging.info("wavelength %f, index %f" %
             (wave_blue, water.getIndex(None, wave_blue).real))

elem.addMaterial("water", water)

elem.addSurface("stop", stopsurf, (None, None))
elem.addSurface("surf1", frontsurf, (None, "water"))
elem.addSurface("surf2", rearsurf, ("water", "water"))
elem.addSurface("surf4", frontsurf, ("water", None))
elem.addSurface("image", image, (None, None))

s.addElement("droplet", elem)

sysseq = [("droplet",
           [
                ("stop", {"is_stop": True}),
                ("surf1", {}),
                ("surf2", {"is_mirror": True}),
                ("surf4", {}),
                ("image", {})])]

raysdict = {
     "radius": dropletradius*0.05,
     "starty": dropletradius*0.9,
     "anglex": -12.*degree,
     "raster": raster.MeridionalFan()
     }

r_red = raytrace(s, sysseq, 11, raysdict, wave=wave_red)[0]
r_blue = raytrace(s, sysseq, 11, raysdict, wave=wave_blue)[0]

draw(s, [(r_red, "red"), (r_blue, "blue")])

system_dump = Serializer(s).serialization
system_gui_toplevel = UIInterfaceClassWithOptimizableVariables(
        s.elements["droplet"].surfaces["surf4"].shape).query_for_dictionary()

#pprint(system_gui_toplevel)
#pprint(system_dump)

fp = open("rainbow.yaml", "wt")
yaml.dump(system_dump, fp)
fp.close()


fp = open("rainbow.json", "wt")
json.dump(system_dump, fp, indent=4)
fp.close()


