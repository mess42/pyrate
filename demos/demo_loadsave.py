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

import numpy as np

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

from pyrateoptics.raytracer.surface_shape import (Asphere,
                                                  Conic,
                                                  ZernikeFringe,
                                                  GridSag)
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

def create_save_load_compare(name, mycreate, mycompare):
    o = mycreate()
    o_dump = Serializer(o).serialization
    save_yaml(mypath + name + ".yaml", o_dump)
    save_json(mypath + name + ".json", o_dump)
    del o_dump
    o_dump2 = load_json(mypath + name + ".json")
    o2 = Deserializer(o_dump2, True, True).class_instance
    assert mycompare(o, o2)


def compare_surface_shapes(ssh1, ssh2):
    x = np.random.random(100)
    y = np.random.random(100)

    z1 = ssh1.getSag(x, y)
    z2 = ssh2.getSag(x, y)

    # TODO: compare also derivatives

    return np.allclose(z1, z2)


def loadsave_conic():
    def create():
        lc = LocalCoordinates.p(name="global")
        c = Conic.p(lc, curv=0.01, cc=-0.1)
        return c

    create_save_load_compare("conic", create,
                             compare_surface_shapes)


def loadsave_asphere():
    def create():
        lc = LocalCoordinates.p(name="global")
        a = Asphere.p(lc, curv=0.01, cc=-1.0, coefficients=[1., 2., 3.])
        return a

    create_save_load_compare("asphere", create,
                             compare_surface_shapes)


def loadsave_gridsag():
    def create():
        lc = LocalCoordinates.p(name="global")
        x = np.linspace(-1, 1, 100)
        (X, Y) = np.meshgrid(x, x)
        Z = X**2 + Y**2
        g = GridSag.p(lc, (x, x, Z))
        return g

    create_save_load_compare("gridsag", create,
                             compare_surface_shapes)


def loadsave_zernike_fringe():
    pass


def loadsave_zernike_standard():
    pass


def loadsave_surface():
    pass


def loadsave_opticalsystem_empty():

    def create():
        s = OpticalSystem.p()
        return s

    def compare(s1, s2):
        # TODO: compare drawing?
        return True

    create_save_load_compare("system_empty", create, compare)


if __name__ == "__main__":
    try:
        os.mkdir(mypath)
    except FileExistsError:
        pass
    # surface shapes
    loadsave_conic()
    loadsave_asphere()
    loadsave_gridsag()
    loadsave_zernike_fringe()
    loadsave_zernike_standard()
    # surface
    loadsave_surface()
    # optical system
    loadsave_opticalsystem_empty()
    # materials
    # isotropic, catalog, anisotropic, grin
