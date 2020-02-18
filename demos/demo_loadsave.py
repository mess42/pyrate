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
    o_dump2_json = load_json(mypath + name + ".json")
    o_dump2_yaml = load_yaml(mypath + name + ".yaml")
    o2_json = Deserializer(o_dump2_json, True, True).class_instance
    o2_yaml = Deserializer(o_dump2_yaml, True, True).class_instance
    return mycompare(o, o2_json) and mycompare(o, o2_yaml)


def compare_surface_shapes(ssh1, ssh2):
    x = np.random.random(100)
    y = np.random.random(100)

    z1 = ssh1.getSag(x, y)
    z2 = ssh2.getSag(x, y)

    # TODO: compare also derivatives

    return np.allclose(z1, z2)


def loadsave_conic():
    def create():
        from pyrateoptics.raytracer.surface_shape import Conic
        lc = LocalCoordinates.p(name="global")
        c = Conic.p(lc, curv=0.01, cc=-0.1)
        return c

    return create_save_load_compare("conic", create,
                                    compare_surface_shapes)


def loadsave_asphere():
    def create():
        from pyrateoptics.raytracer.surface_shape import Asphere
        lc = LocalCoordinates.p(name="global")
        a = Asphere.p(lc, curv=0.01, cc=-1.0, coefficients=[1., 2., 3.])
        return a

    return create_save_load_compare("asphere", create,
                                    compare_surface_shapes)


def loadsave_gridsag():
    def create():
        from pyrateoptics.raytracer.surface_shape import GridSag
        lc = LocalCoordinates.p(name="global")
        x = np.linspace(-1, 1, 100)
        (X, Y) = np.meshgrid(x, x)
        Z = X**2 + Y**2
        g = GridSag.p(lc, (x, x, Z))
        return g

    return create_save_load_compare("gridsag", create,
                                    compare_surface_shapes)


def loadsave_zernike_fringe():
    def create():
        from pyrateoptics.raytracer.surface_shape import ZernikeFringe
        lc = LocalCoordinates.p(name="global")
        x = np.linspace(-1, 1, 100)
        (X, Y) = np.meshgrid(x, x)
        Z = ZernikeFringe.p(lc, normradius=1.0,
                            coefficients=[1., 2., 3., 4., 5., 6., 7., 8.])
        return Z

    return create_save_load_compare("zernikefringe", create,
                                    compare_surface_shapes)


def loadsave_zernike_ansi():
    def create():
        from pyrateoptics.raytracer.surface_shape import ZernikeANSI
        lc = LocalCoordinates.p(name="global")
        x = np.linspace(-1, 1, 100)
        (X, Y) = np.meshgrid(x, x)
        Z = ZernikeANSI.p(lc, normradius=1.0,
                          coefficients=[1., 2., 3., 4., 5., 6., 7., 8.])
        return Z

    return create_save_load_compare("zernikeansi", create,
                                    compare_surface_shapes)


def loadsave_surface():
    def create():
        from pyrateoptics.raytracer.surface_shape import Conic
        from pyrateoptics.raytracer.surface import Surface
        from pyrateoptics.raytracer.aperture import CircularAperture

        lc = LocalCoordinates.p(name="----LCglobal----")
        lc2 = LocalCoordinates.p(name="----LCap----", tiltx=0.1, decx=0.2)
        lc3 = LocalCoordinates.p(name="----LCsh----", tiltx=-0.1, decx=-0.2)

        lc.addChild(lc2)
        lc.addChild(lc3)

        # TODO: child coordinate systems lead to core dumps

        ap = CircularAperture(lc2)
        sh = Conic.p(lc3, curv=0.01, cc=-1)

        su = Surface.p(lc, sh, ap, name="mysurface")

        return su

    def compare(su1, su2):
        sh1 = su1.shape
        sh2 = su2.shape

        cmp_ssh = compare_surface_shapes(sh1, sh2)
        return cmp_ssh

    return create_save_load_compare("surface", create, compare)


def loadsave_opticalsystem_empty():

    def create():
        s = OpticalSystem.p()
        return s

    def compare(s1, s2):
        # TODO: compare drawing?
        return True

    return create_save_load_compare("system_empty", create, compare)


if __name__ == "__main__":

    try:
        os.mkdir(mypath)
    except FileExistsError:
        pass
    my_func_list = [
        # surface shapes
        "loadsave_conic()",
        "loadsave_asphere()",
        "loadsave_gridsag()",
        "loadsave_zernike_fringe()",
        "loadsave_zernike_ansi()",
        # surface
        "loadsave_surface()",
        # optical system
        "loadsave_opticalsystem_empty()"
        # materials
        # isotropic, catalog, anisotropic, grin
    ]

    true_false_list = [eval(fun) for fun in my_func_list]

    for funname, funval in zip(my_func_list, true_false_list):
        print(funname + ": " + str(funval))
