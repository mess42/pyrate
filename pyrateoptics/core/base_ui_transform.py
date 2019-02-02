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
from ..raytracer.globalconstants import degree


def deg2rad(x):
    return x*degree


def rad2deg(x):
    return x/degree


def curv2radius(x):
    return 1./x


def radius2curv(x):
    mycurv = 0.
    if abs(x) > 1e-16:
        mycurv = 1./x
    return mycurv


"""
Provides transformation dict in the form
{"kind1": {"var_name": ("shown_varname", transform_func,
                        inverse_transform_func), ...}, ...}
where transform_func(x) transforms the value into the written
value into UI and inverse_transform_func(x) transforms the value back
e.g.: {"shape_conic": {"curv": ("radius", lambda x: 1./x, lambda x: 1./x)}}
or: {"localcoordinates": {"tiltx": ("tiltx_deg", lambda x: x/degree,
                                    lambda x: x*degree)}}
"""
transformation_dictionary =\
    {"shape_conic": {"curv": ("radius", curv2radius, radius2curv)},
     "localcoordinates": {"tiltx": ("tiltx_deg", rad2deg, deg2rad),
                          "tilty": ("tilty_deg", rad2deg, deg2rad),
                          "tiltz": ("tiltz_deg", rad2deg, deg2rad)}}
