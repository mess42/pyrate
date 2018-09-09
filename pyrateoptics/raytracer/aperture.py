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

import numpy as np

from ..core.log import BaseLogger


class BaseAperture(BaseLogger):
    # for optimizable aperture it would be good to
    # derive from optimizable class, but I think
    # apertures will be inserted after optimization
    """
    Base class representing the aperture of a surface.
    Subclasses may define the actual shapes (circular,
    elliptic, rectangular, etc.)

    Decentering of apertures do not make any sense, since
    the local coordinate system can be redefined if desired.

    The base class does not limit the beam diameter.
    """
    def __init__(self, lc, **kwargs):
        super(BaseAperture, self).__init__(**kwargs)
        self.lc = lc
        self.typicaldimension = 1e16

    def getTypicalDimension(self):
        """
        Returns typical dimension of aperture

        :return self.typicaldimension: float
        """

        return self.typicaldimension

    def arePointsInAperture(self, x, y):
        """

        Returns of points given by numpy arrays x, y are within aperture

        :param x: x in local coordinate system (1xn numpy array of floats)
        :param y: y in local coordinate system (1xn numpy array of floats)

        :return True (1d numpy array of n bools)
        """

        bool_func = self.getBooleanFunction()
        return bool_func(x, y)  # return true always

    def getBooleanFunction(self):
        """

        Returns boolean function of aperture

        :return anonymous function for 2 arguments x, y

        """

        return (lambda x, y: np.ones_like(x, dtype=bool))


class CircularAperture(BaseAperture):
    """
    Circular aperture of a surface.
    """

    def __init__(self, lc, maxradius=1.0, minradius=0.0, **kwargs):
        super(CircularAperture, self).__init__(lc, **kwargs)
        self.maxradius = maxradius
        self.minradius = minradius
        self.typicaldimension = self.maxradius

    def getBooleanFunction(self):
        return (lambda x, y: (x**2 + y**2 >= self.minradius**2) *
                             (x**2 + y**2 <= self.maxradius**2))

    def getDictionary(self):
        res = super(CircularAperture, self).getDictionary()
        res["type"] = "CircularAperture"
        res["minradius"] = self.minradius
        res["maxradius"] = self.maxradius
        return res


class RectangularAperture(BaseAperture):
    """
    Rectangular aperture of a surface.
    """

    def __init__(self, lc, width=1.0, height=1.0, **kwargs):
        super(RectangularAperture, self).__init__(lc, **kwargs)
        self.width = width
        self.height = height
        self.typicaldimension = np.sqrt(self.width**2 + self.height**2)

    def getBooleanFunction(self):
        return (lambda x, y: (x >= -self.width * 0.5) *
                (x <= self.width * 0.5) *
                (y >= -self.height * 0.5) *
                (y <= self.height * 0.5))

    def getDictionary(self):
        res = super(CircularAperture, self).getDictionary()
        res["type"] = "RectangularAperture"
        res["width"] = self.width
        res["height"] = self.height
        return res


accessible_apertures = {None: BaseAperture,
                        "CircularAperture": CircularAperture,
                        "RectangularAperture": RectangularAperture}


def createAperture(lc, ap_dict):

    ap_type = ap_dict.pop("type", None)
    return accessible_apertures[ap_type](lc, **ap_dict)
