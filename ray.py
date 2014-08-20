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
        self.t = zeros(shape(o)[1]) # Geometrical path length to the ray final position.
        self.wave = wave
        self.pol = pol

    def setRayDir(self,k):
        """
        Calculates the unit direction vector of a ray from its wavevector.
        """
        rayDir = 1. * k # copy k, dont just create a pointer
        absk = sqrt( sum(rayDir**2, axis=0) )
        rayDir[0] = rayDir[0] / absk
        rayDir[1] = rayDir[1] / absk
        rayDir[2] = rayDir[2] / absk
        self.rayDir = rayDir

    def draw2d(self, ax, offset=[0,0], color="blue"):
        nrays = shape(self.o)[1]
        for i in arange(nrays):
            y = array( [self.o[1,i], self.o[1,i] + self.t[i] * self.rayDir[1,i] ] )
            z = array( [self.o[2,i], self.o[2,i] + self.t[i] * self.rayDir[2,i] ] )     
            ax.plot(z+offset[1],y+offset[0], color)

 
class RayPath(object):
    def __init__(self, initialraybundle, opticalSystem):
        """
        Class representing the Path of a RayBundle through the whole optical system.

        :param initialraybundle: Raybundle at initial position in the optical system ( RayBundle object )
        :param opticalSystem:  optical system through which the rays are propagated ( OpticalSystem object )

        """
        self.raybundles = [ initialraybundle ]
        N = opticalSystem.getNumberOfSurfaces()

        for i in arange(N-1)+1:
            self.traceToNextSurface(opticalSystem.surfaces[i], opticalSystem.surfaces[i-1].getThickness() )

    def traceToNextSurface(self, nextSurface, thicknessOfCurrentSurface):
        """
        Private routine that propagates a ray bundle to the next surface.
        Please respect the privacy of this class and call it only from methods inside this class.

        :param nextSurface: (Surface object)
        :param thicknessOfCurrentSurface: on-axis geometrical distance between current and nextSurface (float)
        """
        self.raybundles[-1].o[2] -= thicknessOfCurrentSurface 
        intersection, t, normal, validIndices = nextSurface.shap.intersect(self.raybundles[-1])
        self.raybundles[-1].t = t       
        self.raybundles.append(  nextSurface.mater.refract( self.raybundles[-1], intersection, normal, validIndices)  )

    def draw2d(self, opticalsystem, ax, offset=[0,0], color="blue"):
        """
        Plots the surface in a matplotlib figure.
        :param ax: matplotlib subplot handle 
        :param offset: y and z offset (list or 1d numpy array of 2 floats)
        :param vertices: number of points the polygon representation of the surface contains (int)
        :param color: surface draw color (str)
        """
        Nsurf = len(self.raybundles)
        offy = offset[0]
        offz = offset[1]
        for i in arange(Nsurf):
            offz += opticalsystem.surfaces[i].getThickness()
            self.raybundles[i].draw2d(ax, offset = [offy, offz], color = color)
