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


import pickle, pickletools

from core import material
from core import surfShape
from core.optical_system import OpticalSystem, Surface

from core.aperture import CircularAperture, BaseAperture
from core.coordinates import LocalCoordinates

import math

# definition of optical system
s = OpticalSystem() 

lc1 = s.addLocalCoordinateSystem(LocalCoordinates(name="surf1", decz=2.0)) # objectDist
lc2 = s.addLocalCoordinateSystem(LocalCoordinates(name="surf2", decz=3.0))
lc3 = s.addLocalCoordinateSystem(LocalCoordinates(name="surf3", decz=5.0, tiltx=2.5*math.pi/180.0))
lc4 = s.addLocalCoordinateSystem(LocalCoordinates(name="surf4", decz=3.0))
lc5 = s.addLocalCoordinateSystem(LocalCoordinates(name="surf5", decz=3.0))
lc6 = s.addLocalCoordinateSystem(LocalCoordinates(name="surf6", decz=2.0))
lc7 = s.addLocalCoordinateSystem(LocalCoordinates(name="surf7", decz=3.0))
lc8 = s.addLocalCoordinateSystem(LocalCoordinates(name="image", decz=19.0))


s.insertSurface(1, Surface(lc1, surfShape.Conic(curv=1/-5.922), # thickness=3.0,
                           material=material.ConstantIndexGlass(1.7), 
                            aperture=BaseAperture()))

s.insertSurface(2, Surface(lc2, surfShape.Conic(curv=1/-3.160), # thickness=5.0, 
                           aperture=BaseAperture()))

s.insertSurface(3, Surface(lc3, surfShape.Conic(curv=1/15.884), #thickness=3.0,
                           material=material.ConstantIndexGlass(1.7), 
                            aperture=BaseAperture()))

s.insertSurface(4, Surface(lc4, surfShape.Conic(curv=1/-12.756), #thickness=3.0,
                           aperture=BaseAperture()))

s.insertSurface(5, Surface(lc5, surfShape.Conic(), #thickness=2.0, 
                           aperture=BaseAperture())) # Stop Surface

s.insertSurface(6, Surface(lc6, surfShape.Conic(curv=1/3.125), #thickness=3.0,
                           material=material.ConstantIndexGlass(1.5), 
                            aperture=BaseAperture()))

s.insertSurface(7, Surface(lc7, surfShape.Conic(curv=0.1*1/1.479), #thickness=19.0,
                           aperture=BaseAperture()))


s.insertSurface(8, Surface(lc8)) # image

# reintroduce circular apertures

s.surfaces[1].aperture = CircularAperture(0.55)
s.surfaces[2].aperture = CircularAperture(1.0)
s.surfaces[3].aperture = CircularAperture(1.3)
s.surfaces[4].aperture = CircularAperture(1.3)
s.surfaces[5].aperture = CircularAperture(1.01)
s.surfaces[6].aperture = CircularAperture(1.0)
s.surfaces[7].aperture = CircularAperture(1.0)


print "pickle dump"
picklestr = pickle.dumps(s)
pickletools.dis(picklestr)


with open('optical_sys.pkl', 'wb') as output:
    pickle.dump(s, output, pickle.HIGHEST_PROTOCOL)

# print introspection variables
for (s_dict_keys, s_dict_items) in s.__dict__.iteritems():
    print s_dict_keys, " ", s_dict_items

# WARNING: code is operational, but not tested
