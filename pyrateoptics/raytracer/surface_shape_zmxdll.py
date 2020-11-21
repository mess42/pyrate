#!/usr/bin/env python3
# -*- coding: utf-8 -*-
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

import ctypes

from ctypes import c_double, c_byte, c_int, c_void_p, c_char_p

import numpy as np

from ..core.optimizable_variable import FloatOptimizableVariable, FixedState
from .surface_shape import Conic


class USER_DATA(ctypes.Structure):
    _fields_ = [
        ("x", c_double), ("y", c_double), ("z", c_double),
        ("l", c_double), ("m", c_double), ("n", c_double),
        ("ln", c_double), ("mn", c_double), ("nn", c_double),
        ("path", c_double), ("sag1", c_double), ("sag2", c_double),
        ("index", c_double),
        ("dndx", c_double), ("dndy", c_double), ("dndz", c_double),
        ("rel_surf_tran", c_double),
        ("udreserved1", c_double),
        ("udreserved2", c_double),
        ("udreserved3", c_double),
        ("udreserved4", c_double),
        ("string", 20*c_byte)]


class FIXED_DATA(ctypes.Structure):
    _fields_ = [
            ("type", c_int),
            ("numb", c_int),
            ("surf", c_int),
            ("wave", c_int),
            ("wavelength", c_double),
            ("pwavelength", c_double),
            ("n1", c_double),
            ("n2", c_double),
            ("cv", c_double),
            ("thic", c_double),
            ("sdia", c_double),
            ("k", c_double),
            ("param", 9*c_double),
            ("fdreserved1", c_double),
            ("xdata", 201*c_double),
            ("glass", 21*c_byte)
        ]


class FIXED_DATA2(ctypes.Structure):
    _fields_ = [
            ("type", c_int),
            ("numb", c_int),
            ("surf", c_int),
            ("wave", c_int),
            ("unit", c_int),
            ("wavelength", c_double),
            ("pwavelength", c_double),
            ("n1", c_double),
            ("n2", c_double),
            ("cv", c_double),
            ("thic", c_double),
            ("sdia", c_double),
            ("k", c_double),
            ("ireserved", 20*c_int),
            ("dreserved", 20*c_double),
            ("param", 51*c_double),
            ("xdata", 201*c_double),
            ("glass", 21*c_byte)
        ]


class FIXED_DATA3(ctypes.Structure):
    _fields_ = [
            ("type", c_int),
            ("numb", c_int),
            ("surf", c_int),
            ("wave", c_int),
            ("unit", c_int),
            ("serial", c_int),
            ("is_a_mirror", c_int),
            ("is_mirror_space", c_int),
            ("is_air", c_int),
            ("ireserved", 100*c_int),
            ("did_polar", c_int),
            ("max_parameter", c_int),
            ("max_extradata", c_int),
            ("Exr", c_double),
            ("Exi", c_double),
            ("Eyr", c_double),
            ("Eyi", c_double),
            ("Ezr", c_double),
            ("Ezi", c_double),
            ("Ewr", c_double),
            ("Ewi", c_double),
            ("dbreserved", 100*c_double),
            ("wavelength", c_double),
            ("pwavelength", c_double),
            ("n1", c_double),
            ("n2", c_double),
            ("cv", c_double),
            ("thic", c_double),
            ("sdia", c_double),
            ("k", c_double),
            ("param", 201*c_double),
            ("xdata", 501*c_double),
            ("glass", 200*c_byte),
            ("comment", 200*c_byte),
            ("int_data", c_void_p),
            ("db_data", c_void_p),
            ("c_data", c_char_p)
        ]


class FIXED_DATA4(ctypes.Structure):
    _fields_ = [
            ("type", c_int),
            ("numb", c_int),
            ("surf", c_int),
            ("wave", c_int),
            ("unit", c_int),
            ("serial", c_int),
            ("is_a_mirror", c_int),
            ("is_mirror_space", c_int),
            ("is_air", c_int),
            ("ireserved", 100*c_int),
            ("did_polar", c_int),
            ("max_parameter", c_int),
            ("max_extradata", c_int),
            ("Exr", c_double),
            ("Exi", c_double),
            ("Eyr", c_double),
            ("Eyi", c_double),
            ("Ezr", c_double),
            ("Ezi", c_double),
            ("Ewr", c_double),
            ("Ewi", c_double),
            ("dbreserved", 100*c_double),
            ("wavelength", c_double),
            ("pwavelength", c_double),
            ("n1", c_double),
            ("n2", c_double),
            ("cv", c_double),
            ("thic", c_double),
            ("sdia", c_double),
            ("k", c_double),
            ("param", 201*c_double),
            ("xdata", 501*c_double),
            ("glass", 200*c_byte),
            ("comment", 200*c_byte),
            ("data_surf_numvalues", c_int),
            ("db_data", 4801*c_double)
        ]


class FIXED_DATA5(ctypes.Structure):
    _fields_ = [
            ("type", c_int),
            ("numb", c_int),
            ("surf", c_int),
            ("wave", c_int),
            ("unit", c_int),
            ("serial", c_int),
            ("is_a_mirror", c_int),
            ("is_mirror_space", c_int),
            ("is_air", c_int),
            ("ireserved", 100*c_int),
            ("did_polar", c_int),
            ("max_parameter", c_int),
            ("max_extradata", c_int),
            ("Exr", c_double),
            ("Exi", c_double),
            ("Eyr", c_double),
            ("Eyi", c_double),
            ("Ezr", c_double),
            ("Ezi", c_double),
            ("Ewr", c_double),
            ("Ewi", c_double),
            ("dbreserved", 100*c_double),
            ("wavelength", c_double),
            ("pwavelength", c_double),
            ("n1", c_double),
            ("n2", c_double),
            ("cv", c_double),
            ("thic", c_double),
            ("sdia", c_double),
            ("k", c_double),
            ("param", 256*c_double),
            ("glass", 200*c_byte),
            ("comment", 200*c_byte),
            ("data_surf_numvalues", c_int),
            ("db_data", 4801*c_double)
        ]


fixed_data_dict = {"": FIXED_DATA,
                   "2": FIXED_DATA2,
                   "3": FIXED_DATA3,
                   "4": FIXED_DATA4,
                   "5": FIXED_DATA5}


class ZMXDLLShape(Conic):
    """
    This surface is able to calculate certain
    shape quantities from a DLL loaded externally.
    """

    @classmethod
    def p(cls, lc, dllfile,
          param_dict=None,
          xdata_dict=None,
          is_win_dll=False,
          max_uds_index_check=10,
          curv=0.0,
          cc=0.0, name=""):

        """
        param: lc LocalCoordinateSystem of shape
        param: dllfile (string) path to DLL file
        param: param_dict (dictionary) key: string, value:
            (int 0 to 8, float); initializes optimizable variables.
        param: xdata_dict (dictionary) key: string, value:
            (int  0 to 200, float); initializes optimizable variables.
        param: is_win_dll (bool): is DLL compiled in Windows or Linux?

        Notice that this class only calls the appropriate
        functions from the DLL. The user is responsible to
        get all necessary functions running to use this DLL.
        (i.e. intersect, sag, normal, ...)

        To compile the DLL for Linux, remove all Windows and
        calling convention stuff and build it
        via:

            gcc -c -fpic -o us_stand.o us_stand.c -lm
            gcc -shared -o us_stand.so us_stand.o

        via MingW it is possible to compile for Windows:

            For a 32bit DLL:
            i686-w64-mingw32-gcc -shared -o us_stand.dll us_stand.c

            For a 64bit DLL:
            x86_64-w64-mingw32-gcc -shared -o us_stand.dll us_stand.c


        """

        if param_dict is None:
            param_dict = {}

        if xdata_dict is None:
            xdata_dict = {}

        zmx_dll_structure = {}
        zmx_dll_annotations = {}

        zmx_dll_structure["lc"] = lc
        zmx_dll_structure["curvature"] = FloatOptimizableVariable(
            FixedState(curv), name="curv")
        zmx_dll_structure["cc"] = FloatOptimizableVariable(
            FixedState(cc), name="conic constant")

        if is_win_dll:
            mydll = ctypes.WinDLL(dllfile)
        else:
            mydll = ctypes.CDLL(dllfile)

        zmx_dll_structure["param"] = {}
        for (key, (value_int, value_float)) in param_dict.items():
            zmx_dll_structure["param"][value_int] = FloatOptimizableVariable(
                FixedState(value_float),
                name="param" + str(value_int))
        zmx_dll_structure["xdata"] = {}
        for (key, (value_int, value_float)) in xdata_dict.items():
            zmx_dll_structure["xdata"][value_int] = FloatOptimizableVariable(
                FixedState(value_float),
                name="xdata" + str(value_int))
        # check for UserDefinedSurfaceX where X is the index
        # For the time being, this wrapper class supports only
        # one surface function (either UserDefinedSurface or
        # UserDefinedSurface2 or ...)
        final_uds_index = None
        for uds_index in range(max_uds_index_check):
            try:
                check_uds_index = str(uds_index) if uds_index > 0 else ""
                us_surf = getattr(mydll, "UserDefinedSurface" +
                                  check_uds_index)
            except AttributeError:
                pass
            else:
                final_uds_index = check_uds_index

        zmx_dll_structure["dll"] = mydll
        zmx_dll_structure["us_surf"] =\
            us_surf
        zmx_dll_annotations["us_index"] = final_uds_index

        return cls(zmx_dll_annotations, zmx_dll_structure, name=name)

    def setKind(self):
        self.kind = "shape_ZMXDLLShape"

    def write_param(self, fixed_data):
        for (key, var) in self.param.items():
            fixed_data.param[key] = var()
        return fixed_data

    def write_xdata(self, fixed_data):
        for (key, var) in self.xdata.items():
            fixed_data.xdata[key] = var()
        return fixed_data

    def intersect(self, raybundle):
        (r0, ray_dir) = self.getLocalRayBundleForIntersect(raybundle)

        intersection = np.zeros_like(r0)

        user_data = USER_DATA()
        fixed_data = fixed_data_dict[self.annotations["us_index"]]()

        fixed_data.type = 5 # ask for intersection
        fixed_data.k = self.cc()
        fixed_data.cv = self.curvature()
        fixed_data.wavelength = raybundle.wave

        fixed_data = self.write_param(fixed_data)
        fixed_data = self.write_xdata(fixed_data)

        (x0, y0, z0) = r0.T

        (l_cos, m_cos, n_cos) = ray_dir.T

        myvalues = np.vstack((x0, y0, z0, l_cos, m_cos, n_cos)).T

        for (ind, (xl, yl, zl, ll, ml, nl)) in enumerate(myvalues.tolist()):

            user_data.x = xl
            user_data.y = yl
            user_data.z = zl
            user_data.l = ll
            user_data.m = ml
            user_data.n = nl

            self.us_surf(ctypes.byref(user_data),
                         ctypes.byref(fixed_data))

            intersection[0, ind] = user_data.x
            intersection[1, ind] = user_data.y
            intersection[2, ind] = user_data.z

        globalinter = self.lc.returnLocalToGlobalPoints(intersection)

        raybundle.append(globalinter, raybundle.k[-1], raybundle.Efield[-1],
                         raybundle.valid[-1])

    def getSag(self, x, y):
        user_data = USER_DATA()
        fixed_data = fixed_data_dict[self.annotations["us_index"]]()

        fixed_data.type = 3  # ask for sag
        fixed_data.k = self.cc()
        fixed_data.cv = self.curvature()

        fixed_data = self.write_param(fixed_data)
        fixed_data = self.write_xdata(fixed_data)

        z = np.zeros_like(x)
        myvalues = np.vstack((x, y)).T

        for (ind, (xp, yp)) in enumerate(myvalues.tolist()):
            user_data.x = xp
            user_data.y = yp
            retval = self.us_surf(ctypes.byref(user_data),
                                  ctypes.byref(fixed_data))
            z[ind] = user_data.sag1 if retval == 0 else user_data.sag2

        return z
