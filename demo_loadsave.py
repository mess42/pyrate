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

from core.aperture import CircularAperture

# definition of optical system
# definition of optical system
# definition of optical system
s = OpticalSystem(objectDistance = 2.0)

s.insertSurface(1, Surface(surfShape.Conic(curv=1/-5.922), thickness=3.0,
                           material=material.ConstantIndexGlass(1.7), aperture=CircularAperture(0.55)))
s.insertSurface(2, Surface(surfShape.Conic(curv=1/-3.160), thickness=5.0, aperture=CircularAperture(1.0)))
s.insertSurface(3, Surface(surfShape.Conic(curv=1/15.884), thickness=3.0,
                           material=material.ConstantIndexGlass(1.7), aperture=CircularAperture(1.3)))
s.insertSurface(4, Surface(surfShape.Conic(curv=1/-12.756), thickness=3.0,
                           aperture=CircularAperture(1.3)))
s.insertSurface(5, Surface(surfShape.Conic(), thickness=2.0, aperture=CircularAperture(1.01))) # Stop Surface
s.insertSurface(6, Surface(surfShape.Conic(curv=1/3.125), thickness=3.0,
                           material=material.ConstantIndexGlass(1.5), aperture=CircularAperture(1.0)))
s.insertSurface(7, Surface(surfShape.Conic(curv=1/1.479), thickness=19.0,
                           aperture=CircularAperture(1.0)))

print "pickle dump"
picklestr = pickle.dumps(s)
pickletools.dis(picklestr)


with open('optical_sys.pkl', 'wb') as output:
    pickle.dump(s, output, pickle.HIGHEST_PROTOCOL)

# print introspection variables
for (s_dict_keys, s_dict_items) in s.__dict__.iteritems():
    print s_dict_keys, " ", s_dict_items

