#!/usr/bin/env/python
"""
Pyrate - Optical raytracing based on Python

Copyright (C) 2014 Moritz Esslinger moritz.esslinger@web.de
               and    Johannes Hartung j.hartung@gmx.net
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
    def __init__(self, tx=0.0, ty=0.0):
        self.typicaldimension = 1e10
        self.tx = tx
        self.ty = ty

    def getTypicalDimension(self):
        return self.typicaldimension

    def arePointsInAperture(self, x, y):
        return np.ones_like(x, dtype=bool) # return true always

    def getBooleanFunction(self):
        return (lambda x, y: True)


class CircularAperture(BaseAperture):
    """
    Circular aperture of a surface.
    """

    def __init__(self, semidiameter = 1.0, tx = 0.0, ty = 0.0):
        super(CircularAperture, self).__init__(tx, ty)
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

    def __init__(self, w=1.0, h=1.0, tx=0.0, ty=0.0):
        super(RectangularAperture, self).__init__(tx, ty)
        self.width = w
        self.height = h
        self.typicaldimension = math.sqrt(self.width**2 + self.height**2)



    def arePointsInAperture(self, x, y):
        return x >= -self.width*0.5 - self.tx and x <= self.width*0.5 - self.tx and \
            y >= -self.height*0.5 - self.ty and y <= self.height*0.5 - self.ty

    def getBooleanFunction(self):
        return (lambda x, y: x >= -self.width*0.5 - self.tx and x <= self.width*0.5 - self.tx and \
            y >= -self.height*0.5 - self.ty and y <= self.height*0.5 - self.ty)



