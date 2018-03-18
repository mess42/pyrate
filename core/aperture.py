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
import math

# TODO: method headers

class BaseAperture(object):
    # for optimizable aperture it would be good to derive from optimizable class, but I think
    # apertures will be inserted after optimization
    """
    Base class representing the aperture of a surface.
    Subclasses may define the actual shapes (circular,
    elliptic, rectangular, etc.)

    The base class does not limit the beam diameter.
    """
    def __init__(self, lc, tx=0.0, ty=0.0):
        self.lc = lc
        self.typicaldimension = 1e10
        self.tx = tx
        self.ty = ty

    def getTypicalDimension(self):
        """
        Returns typical dimension of aperture

        :return self.typicaldimension: float
        """

        return self.typicaldimension

    def arePointsInAperture(self, x, y):
        """

        Returns of points given by numpy arrays x, y are within aperture

        :param x: x location of points in local coordinate system (1d numpy array of n floats)
        :param y: y location of points in local coordinate system (1d numpy array of n floats)

        :return True (1d numpy array of n bools)
        """


        return np.ones_like(x, dtype=bool) # return true always

    def getBooleanFunction(self):
        """

        Returns boolean function of aperture

        :return anonymous function for 2 arguments x, y

        """


        return (lambda x, y: True)


class CircularAperture(BaseAperture):
    """
    Circular aperture of a surface.
    """

    def __init__(self, lc, semidiameter = 1.0, tx = 0.0, ty = 0.0):
        super(CircularAperture, self).__init__(lc, tx, ty)
        self.semidiameter = semidiameter
        self.typicaldimension = self.semidiameter

    def arePointsInAperture(self, x, y):
        return (x - self.tx)**2 + (y - self.ty)**2 <= self.semidiameter**2

    def getBooleanFunction(self):
        return (lambda x, y: (x - self.tx)**2 + (y - self.ty)**2 <= self.semidiameter**2)




class RectangularAperture(BaseAperture):
    """
    Rectangular aperture of a surface.
    """

    def __init__(self, lc, w=1.0, h=1.0, tx=0.0, ty=0.0):
        super(RectangularAperture, self).__init__(lc, tx, ty)
        self.width = w
        self.height = h
        self.typicaldimension = math.sqrt(self.width**2 + self.height**2)



    def arePointsInAperture(self, x, y):
        return (x >= -self.width*0.5 - self.tx)*(x <= self.width*0.5 - self.tx)* \
            (y >= -self.height*0.5 - self.ty)*(y <= self.height*0.5 - self.ty)

    def getBooleanFunction(self):
        return (lambda x, y: (x >= -self.width*0.5 - self.tx)*(x <= self.width*0.5 - self.tx)* \
            (y >= -self.height*0.5 - self.ty)*(y <= self.height*0.5 - self.ty))



