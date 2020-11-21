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

from ..raytracer.globalconstants import degree


def deg2rad(deg):
    """
    Transforms degrees in radians for UI.
    """
    return deg*degree


def rad2deg(rad):
    """
    Transforms radians into degrees for UI.
    """
    return rad/degree


def curv2radius(curv):
    """
    Transforms curvature to radius for UI.
    1e-16 is the limit, where the radius is set to zero.
    """
    myradius = 0.
    if abs(curv) > 1e-16:
        myradius = 1./curv
    return myradius


def radius2curv(radius):
    """
    Transforms radius to curvature for UI.
    1e-16 is the limit, where the curvature is set to zero.
    """
    mycurv = 0.
    if abs(radius) > 1e-16:
        mycurv = 1./radius
    return mycurv


def radius_string(rad_string):
    """
    Transforms radius string into a number from UI.
    """
    result = rad_string
    if rad_string.lower() == "0.0" or rad_string.lower() == "0":
        result = "infinity"
    if rad_string.lower() == "inf" or\
       rad_string.lower() == "oo" or\
       rad_string.lower() == "infinity":
        result = "0.0"
    return result

# Provides transformation dict in the form
# {"kind1": {"var_name": ("shown_varname", transform_func,
#                         inverse_transform_func), ...}, ...}
# where transform_func(x) transforms the value into the written
# value into UI and inverse_transform_func(x) transforms the value back
# e.g.: {"shape_conic": {"curv": ("radius", lambda x: 1./x, lambda x: 1./x)}}
# or: {"localcoordinates": {"tiltx": ("tiltx_deg", lambda x: x/degree,
#                                     lambda x: x*degree)}}

# Translates values into human readable form. (e.g. radius better than curv,
# deg better than rad)

TRANSFORMATION_DICTIONARY_TO_UI =\
    {"shape_Conic": {"curvature": ("radius", curv2radius)},
     "shape_Cylinder": {"curvature": ("radius", curv2radius)},
     "shape_Asphere": {"curv": ("radius", curv2radius)},
     "shape_Biconic": {"curvx": ("radiusx", curv2radius),
                       "curvy": ("radiusy", curv2radius)},
     "localcoordinates": {"tiltx": ("tiltx_deg", rad2deg),
                          "tilty": ("tilty_deg", rad2deg),
                          "tiltz": ("tiltz_deg", rad2deg)}}

# Translates values back into optimization friendly form. (e.g. curv better
# than radius, rad better than deg.)


TRANSFORMATION_DICTIONARY_FROM_UI =\
    {"shape_Conic": {"radius": ("curvature", radius2curv)},
     "shape_Cylinder": {"radius": ("curvature", radius2curv)},
     "shape_Asphere": {"radius": ("curv", radius2curv)},
     "shape_Biconic": {"radiusx": ("curvx", radius2curv),
                       "radiusy": ("curvy", radius2curv)},
     "localcoordinates": {"tiltx_deg": ("tiltx", deg2rad),
                          "tilty_deg": ("tilty", deg2rad),
                          "tiltz_deg": ("tiltz", deg2rad)}}

# Translates converted string values into better readable/writable form.

TRANSFORMATION_DICTIONARY_UI_STRING = {"radius": radius_string}
