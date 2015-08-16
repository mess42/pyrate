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

from numpy import *
import aperture
import material


class RayBundle(object):
    def __init__(self, o, k, rayID, wave=0.55, pol=[]):
        """
        Class representing a bundle of rays.

        :param o:     Origin of the rays.   (2d numpy 3xN array of float)
        :param k:     Wavevectors of the rays, normalized by 2pi/lambda.
                      (2d numpy 3xN array of float)
                      For media obeying the Snell law, the length of each vector
                      is the  refractive index of the current medium.
        :param rayID: Set an ID number for each ray in the bundle (1d numpy array of int)
                      (for example, the ray index at surface 0)
        :param wave:  Wavelength of the radiation in micrometers. (float)
        :param pol:   Polarization state of the rays. (2d numpy 2xN array of complex); not implemented yet

        to do: implement polarization
        to do: implement reflection / transmission coefficients
               coherent (Jones formalism) or incoherent (Mueller/Stokes formalism) ?
               the difference is, that Jones allows for description of phases of the coherent, fully polarized beam
               and Mueller allows for describing intensity transmission of partially polarized beams
        """
        self.o = o
        self.k = k
        self.rayID = rayID
        self.setRayDir(k)
        self.t = zeros(shape(o)[1])  # Geometrical path length to the ray final position.
        self.wave = wave
        self.pol = pol

    def setRayDir(self, k):
        """
        Calculates the unit direction vector of a ray from its wavevector.
        """
        rayDir = 1. * k  # copy k, dont just create a pointer
        absk = sqrt(sum(rayDir**2, axis=0))
        rayDir[0] = rayDir[0] / absk
        rayDir[1] = rayDir[1] / absk
        rayDir[2] = rayDir[2] / absk
        self.rayDir = rayDir

    def getCentroidPos(self):
        """
        Returns the arithmetic average position of all rays at the origin of the ray bundle.

        :return centr: centroid position (1d numpy array of 3 floats)
        """
        oneOverN = 1.0 / (shape(self.o)[1])
        xav = sum(self.o[0][0:]) * oneOverN
        yav = sum(self.o[1][0:]) * oneOverN
        zav = sum(self.o[2][0:]) * oneOverN

        return array([xav, yav, zav])

    def getChiefPos(self):
        """
        Returns the chief ray position at the origin of the ray bundle.

        :return chief: chief position (1d numpy array of 3 floats)
        """
        return self.o[:, 0]

    def getRMSspotSize(self, referencePos):
        """
        Returns the root mean square (RMS) deviation of all ray positions
        with respect to a reference position at the origin of the ray bundle.

        :referencePos: (1d numpy array of 3 floats)

        :return rms: RMS spot size (float)
        """
        deltax = self.o[0][1:] - referencePos[0]
        deltay = self.o[1][1:] - referencePos[1]
        deltaz = self.o[2][1:] - referencePos[2]

        N = len(deltax)

        return sqrt((sum(deltax**2) + sum(deltay**2) + sum(deltaz**2)) / (N-1.0))

    def getRMSspotSizeCentroid(self):
        """
        Returns the root mean square (RMS) deviation of all ray positions
        with respect to the centroid at the origin of the ray bundle.

        :return rms: RMS spot size (float)
        """
        centr = self.getCentroidPos()
        return self.getRMSspotSize(centr)

    def getRMSspotSizeChief(self):
        """
        Returns the root mean square (RMS) deviation of all ray positions
        with respect to the chief ray at the origin of the ray bundle.

        :return rms: RMS spot size (float)
        """
        chief = self.getChiefPos()
        return self.getRMSspotSize(chief)

    def getCentroidDirection(self):
        """
        Returns the arithmetic average direction of all rays at the origin of the ray bundle.

        :return centr: centroid unit direction vector (1d numpy array of 3 floats)
        """
        xav = sum(self.rayDir[0][1:])
        yav = sum(self.rayDir[1][1:])
        zav = sum(self.rayDir[2][1:])

        length = sqrt(xav**2 + yav**2 + zav**2)

        return array([xav, yav, zav]) / length

    def getChiefDirection(self):
        """
        Returns the chief ray unit direction vector.

        :return chief: chief unit direction (1d numpy array of 3 floats)
        """
        return self.rayDir[:, 0]

    def getRMSangluarSize(self, refDir):
        """
        Returns the root mean square (RMS) deviation of all ray directions
        with respect to a reference direction.
        The return value is approximated for small angluar deviations from refDir.

        :param refDir: reference direction vector (1d numpy array of 3 floats)
                       Must be normalized to unit length.

        :return rms: RMS angular size in rad (float)
        """

        # sin(angle between rayDir and refDir) = abs( rayDir crossproduct refDir )
        # sin(angle)**2 approx angle**2 for small deviations from the reference,
        # but for large deviations the definition makes no sense, anyway

        crossX = self.rayDir[1][1:] * refDir[2] - self.rayDir[2][1:] * refDir[1]
        crossY = self.rayDir[2][1:] * refDir[0] - self.rayDir[0][1:] * refDir[2]
        crossZ = self.rayDir[0][1:] * refDir[1] - self.rayDir[1][1:] * refDir[0]
        N = len(crossX)

        return sqrt(sum(crossX**2 + crossY**2 + crossZ**2) / (N-1.0))

    def getRMSangluarSizeCentroid(self):
        """
        Returns the root mean square (RMS) deviation of all ray directions
        with respect to the centroid direction.

        :return rms: RMS angular size in rad (float)
        """
        return self.getRMSangluarSize(self.getCentroidDirection())

    def getRMSangluarSizeChief(self):
        """
        Returns the root mean square (RMS) deviation of all ray directions
        with respect to the chief direction.

        :return rms: RMS angular size in rad (float)
        """
        return self.getRMSangluarSize(self.getChiefDirection())

    def draw2d(self, ax, offset=(0, 0), color="blue"):
        nrays = shape(self.o)[1]
        for i in arange(nrays):
            y = array([self.o[1, i], self.o[1, i] + self.t[i] * self.rayDir[1, i]])
            z = array([self.o[2, i], self.o[2, i] + self.t[i] * self.rayDir[2, i]])
            ax.plot(z+offset[1], y+offset[0], color)


class RayPath(object):
    def __init__(self, initialraybundle, opticalSystem):
        """
        Class representing the Path of a RayBundle through the whole optical system.

        :param initialraybundle: Raybundle at initial position in the optical system ( RayBundle object )
        :param opticalSystem:  optical system through which the rays are propagated ( OpticalSystem object )

        """
        self.raybundles = [initialraybundle]
        N = opticalSystem.getNumberOfSurfaces()

        for i in arange(N-1)+1:
            self.traceToNextSurface(opticalSystem.surfaces[i-1], opticalSystem.surfaces[i])

    def traceToNextSurface(self, actualSurface, nextSurface):
        """
        Private routine that propagates a ray bundle to the next surface.
        Thickness can be extracted from actualSurface.
        Please respect the privacy of this class and call it only from methods inside this class.

        :param actualSurface: (Surface object)
        :param nextSurface: (Surface object)
        """
        
        if isinstance(actualSurface.material, material.GrinMaterial):
            intersection, t, normal, validindices, propraybundles = actualSurface.material.propagate(
                                                                                                     nextSurface, 
                                                                                                     self.raybundles[-1]
                                                                                                     )
            self.raybundles += propraybundles
        else:
            # this if-path is for linear ray transfer between some surfaces
            self.raybundles[-1].o[2] -= actualSurface.getThickness()
            intersection, t, normal, validIndices = nextSurface.shape.intersect(self.raybundles[-1])#, aperture.BaseAperture())

            validIndices *= nextSurface.aperture.arePointsInAperture(intersection[0], intersection[1])
            validIndices[0] = True # hail to the chief ray

            # finding valid indices due to an aperture is not in responsibility of the surfShape class anymore
            # TODO: needs heavy testing

            self.raybundles[-1].t = t
        self.raybundles.append(nextSurface.material.refract(self.raybundles[-1], intersection, normal, validIndices))

    def draw2d(self, opticalsystem, ax, offset=(0, 0), color="blue"):
        Nsurf = len(self.raybundles)
        offy = offset[0]
        offz = offset[1]
        for i in arange(Nsurf):
            offz += opticalsystem.surfaces[i].thickness.val
            self.raybundles[i].draw2d(ax, offset=(offy, offz), color=color)

