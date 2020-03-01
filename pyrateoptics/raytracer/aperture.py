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

import math
import numpy as np

from ..core.base import ClassWithOptimizableVariables


class BaseAperture(ClassWithOptimizableVariables):
    # for optimizable aperture it would be good to
    # derive from optimizable class, but I think
    # apertures will be inserted after optimization

    # comment 20200217: class has to be derived from
    # ClassWithOptimizableVariables anyway, for following
    # reasons:
    # * localcoordinate system is still optimizable
    # * it is saved by serializer
    """
    Base class representing the aperture of a surface.
    Subclasses may define the actual shapes (circular,
    elliptic, rectangular, etc.)

    Decentering of apertures do not make any sense, since
    the local coordinate system can be redefined if desired.

    The base class does not limit the beam diameter.
    """
    @classmethod
    def p(cls, lc, name="", *_):
        return cls({"typicaldimension": 1e16}, {"lc": lc},
                   name=name)

    def setKind(self):
        self.kind = "aperture"

    def get_typical_dimension(self):
        """
        Returns typical dimension of aperture

        :return self.typicaldimension: float
        """

        return self.annotations["typicaldimension"]

    def are_points_in_aperture(self, x_intersection, y_intersection):
        """

        Returns of points given by numpy arrays x, y are within aperture

        :param x_intersection: x in local coordinates (1xn npy array of floats)
        :param y_intersection: y in local coordinates (1xn npy array of floats)

        :return True (1d numpy array of n bools)
        """

        bool_func = self.get_boolean_function()
        return bool_func(x_intersection, y_intersection)  # return true always

    def get_boolean_function(self):
        """

        Returns boolean function of aperture

        :return anonymous function for 2 arguments x, y

        """

        return (lambda x, y: np.ones_like(x, dtype=bool))


class CircularAperture(BaseAperture):
    """
    Circular aperture of a surface.
    """

    def setKind(self):
        self.kind = "aperture_Circular"

    @classmethod
    def p(cls, lc, maxradius=1.0, minradius=0.0, name="", *_):
        return cls({"maxradius": maxradius,
                    "minradius": minradius,
                    "typicaldimension": maxradius},
                   {"lc": lc}, name=name)

    def get_boolean_function(self):
        return (lambda x, y:
                (x**2 + y**2 >= self.annotations["minradius"]**2) *
                (x**2 + y**2 <= self.annotations["maxradius"]**2))


class RectangularAperture(BaseAperture):
    """
    Rectangular aperture of a surface.
    """

    def setKind(self):
        self.kind = "aperture_Rectangle"

    @classmethod
    def p(cls, lc, width=1.0, height=1.0, name="", *_):
        cls({"width": width,
             "height": height,
             "typicaldimension": math.sqrt(width**2 + height**2)},
            {"lc": lc}, name=name)

    def get_boolean_function(self):
        width = self.annotations["width"]
        height = self.annotations["height"]
        return (lambda x, y: (x >= -width * 0.5) *
                (x <= width * 0.5) *
                (y >= -height * 0.5) *
                (y <= height * 0.5))


# Needed for convenience functions in pyrateoptics

ACCESSIBLE_APERTURES = {None: BaseAperture,
                        "CircularAperture": CircularAperture,
                        "RectangularAperture": RectangularAperture}


def create_aperture(localcoordinates, ap_dict):
    """
    Creates aperture object from dictionary
    """
    ap_type = ap_dict.pop("type", None)
    return ACCESSIBLE_APERTURES[ap_type].p(localcoordinates, **ap_dict)
