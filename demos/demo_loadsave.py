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

import os
import json
import yaml
import pprint

import logging

# found_jsonpickle = True
# try:
#     import jsonpickle
# except ImportError:
#     found_jsonpickle = False
# print(s)


# if found_jsonpickle:
#     print("pickle dump")
#     frozen = jsonpickle.encode(s)
#
#     with open('optical_sys.jpkl', 'w') as output:
#         output.write(frozen)

# WARNING: code is operational, but not tested

from pyrateoptics.raytracer.surface_shape import Asphere, Conic, ZernikeFringe
from pyrateoptics.raytracer.localcoordinates import LocalCoordinates
from pyrateoptics.raytracer.optical_system import OpticalSystem
from pyrateoptics.core.serializer import Serializer, Deserializer


logging.basicConfig(level=logging.DEBUG)

# definition of optical system

mypath = "serialization_stuff/"


def load_json(filename):
    fp = open(filename, "rt")
    mylist = json.load(fp)
    fp.close()
    return mylist


def save_json(filename, mydump):
    fp = open(filename, "wt")
    json.dump(mydump, fp, indent=4)
    fp.close()


def load_yaml(filename):
    fp = open(filename, "rt")
    mylist = yaml.load(fp)
    fp.close()
    return mylist


def save_yaml(filename, mydump):
    fp = open(filename, "wt")
    yaml.dump(mydump, fp)
    fp.close()


# generalize check to:
#   * create object
#   * save/load object json/yaml
#   * compare reconstructed object with original

def loadsave_conic():
    lc = LocalCoordinates.p(name="global")
    c = Conic.p(lc, curv=0.01, cc=-0.1)
    c_dump = Serializer(c).serialization
    save_yaml(mypath + "conic.yaml", c_dump)
    save_json(mypath + "conic.json", c_dump)

    del c_dump
    c_dump = load_json(mypath + "conic.json")
    c2 = Deserializer(c_dump, True, True).class_instance
    print(c2.getSag(1.0, 0.0) - c.getSag(1., 0.))


def loadsave_asphere():
    lc = LocalCoordinates.p(name="global")
    a = Asphere.p(lc, curv=0.01, cc=-1.0, coefficients=[1., 2., 3.])
    a_dump = Serializer(a).serialization
    save_yaml(mypath + "asphere.yaml", a_dump)
    save_json(mypath + "asphere.json", a_dump)

    del a_dump
    a_dump = load_json(mypath + "asphere.json")
    a2 = Deserializer(a_dump, True, True).class_instance
    print(a2.getSag(1.0, 0.0) - a.getSag(1., 0.))


def loadsave_zernike_fringe():
    pass


def loadsave_zernike_standard():
    pass


def loadsave_surface():
    pass


def loadsave_opticalsystem_empty():

    s = OpticalSystem.p()

    s_dump = Serializer(s).serialization
    save_yaml(mypath + "system_empty.yaml", s_dump)
    save_json(mypath + "system_empty.json", s_dump)

    del s_dump
    s_dump = load_json(mypath + "system_empty.json")
    s2 = Deserializer(s_dump, True, True).class_instance

    print(s)
    print(s2)


if __name__ == "__main__":
    try:
        os.mkdir(mypath)
    except FileExistsError:
        pass
    loadsave_conic()
    loadsave_asphere()
    loadsave_zernike_fringe()
    loadsave_zernike_standard()
    loadsave_surface()
    loadsave_opticalsystem_empty()
    # materials
    # isotropic, catalog, anisotropic, grin
