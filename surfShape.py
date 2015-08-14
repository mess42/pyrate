#!/usr/bin/env/python
"""
Pyrate - Optical raytracing based on Python

Copyright (C) 2014 Moritz Esslinger moritz.esslinger@web.de
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
from optimize import ClassWithOptimizableVariables
from aperture import CircularAperture


class Shape(ClassWithOptimizableVariables):
    def __init__(self):
        """
        Virtual Class for all surface shapes.
        The shape of a surface provides a function to calculate
        the intersection point with a ray.
        """
        super(Shape, self).__init__()

    def intersect(self, raybundle):
        """
        Intersection routine returning intersection point
        with ray and normal vector of the surface.
        :param raybundle: RayBundle that shall intersect the surface. (RayBundle Object)
        :return t: geometrical path length to the next surface (1d numpy array of float)
        :return normal: surface normal vectors (2d numpy 3xN array of float)
        :return validIndices: whether indices hit the surface (1d numpy array of bool)
        """
        raise NotImplementedError()

    def getSag(self, x, y):
        """
        Returns the sag of the surface for given coordinates - mostly used
        for plotting purposes.
        :param x: x coordinate perpendicular to the optical axis (list or numpy 1d array of float)
        :param y: y coordinate perpendicular to the optical axis (list or numpy 1d array of float)
        :return z: sag (list or numpy 1d array of float)
        """
        raise NotImplementedError()

    def getCentralCurvature(self):
        """
        Returns the curvature ( inverse local radius ) on the optical axis.
        :return curv: (float)
        """
        raise NotImplementedError()

    def draw2d(self, ax, offset=(0, 0), vertices=100, color="grey"):
        """
        Plots the surface in a matplotlib figure.
        :param ax: matplotlib subplot handle
        :param offset: y and z offset (list or 1d numpy array of 2 floats)
        :param vertices: number of points the polygon representation of the surface contains (int)
        :param color: surface draw color (str)
        """
        raise NotImplementedError()

    def draw3d(self, offset=(0, 0, 0), tilt=(0, 0, 0), color="grey"):
        """
        To do: find fancy rendering package
        """
        raise NotImplementedError()


class Conic(Shape):
    def __init__(self, curv=0.0, cc=0.0, semidiam=0.0):
        """
        Create rotationally symmetric conic section surface.

        :param curv: Curvature of the surface (float).
        :param cc: Conic constant (float).
        :param semidiam: Semi-diameter of the surface (float).

        -1 < cc < 0 oblate ellipsoid
             cc = 0 sphere
         0 < cc < 1 prolate ellipsoid
             cc = 1 paraboloid
             cc > 1 hyperboloid
        """
        super(Conic, self).__init__()

        self.curvature = self.createOptimizableVariable("curvature", value=curv, status=False)
        self.conic = self.createOptimizableVariable("conic constant", value=cc, status=False)
        #self.sdia = self.createOptimizableVariable("semi diameter", value=semidiam, status=False)

    def getSag(self, x, y):
        rs = x**2 + y**2
        sqrtterm = 1 - (1+self.conic.val) * self.curvature.val**2 * rs
        res =  self.curvature.val * rs / (1 + np.sqrt(sqrtterm))

        return res

    def getCentralCurvature(self):
        # Conic curvature on axis is only influenced by spherical curvature term
        return self.curvature.val

    def intersect(self, raybundle):
        rayDir = raybundle.rayDir

        r0 = raybundle.o

        F = rayDir[2] - self.curvature.val * (rayDir[0] * r0[0] + rayDir[1] * r0[1] + rayDir[2] * r0[2] * (1+self.conic.val))
        G = self.curvature.val * (r0[0]**2 + r0[1]**2 + r0[2]**2 * (1+self.conic.val)) - 2 * r0[2]
        H = - self.curvature.val - self.conic.val * self.curvature.val * rayDir[2]**2

        square = F**2 + H*G

        t = G / (F + np.sqrt(square))

        intersection = r0 + raybundle.rayDir * t

        # find indices of rays that don't intersect with the sphere
        validIndices = (square > 0) #*(intersection[0]**2 + intersection[1]**2 <= 10.0**2))
        # finding valid indices due to an aperture is not in responsibility of the surfShape class anymore
        validIndices[0] = True  # hail to the chief

        # Normal
        normal = np.zeros(np.shape(r0), dtype=float)
        normal[0] = -self.curvature.val * intersection[0]
        normal[1] = -self.curvature.val * intersection[1]
        normal[2] = 1 - self.curvature.val * intersection[2] * (1+self.conic.val)

        absn = np.sqrt(np.sum(normal**2, axis=0))

        normal[0] = normal[0] / absn
        normal[1] = normal[1] / absn
        normal[2] = normal[2] / absn

        return intersection, t, normal, validIndices

    def draw2d(self, ax, offset=(0, 0), vertices=100, color="grey", ap=None):
        if ap == None:
            effdia = 10.0
        else:
            if ap.getTypicalDimension() <= 1000.0:
                # TODO: maybe introduce aperture types Object and Image to distuingish from very large normal apertures
                effdia = ap.getTypicalDimension() #self.sdia.val if self.sdia.val < 10.0 else 10.0
            else:
                effdia = 10.0
        y = effdia * np.linspace(-1, 1, vertices)
        isyap = np.array(ap.arePointsInAperture(np.zeros_like(y), y))
        yinap = y[isyap]
        zinap = self.getSag(0, yinap)
        ax.plot(zinap+offset[1], yinap+offset[0], color)


class Asphere(Shape):
    """
    to do: polynomial asphere as base class for sophisticated surface descriptions
    """
    pass


