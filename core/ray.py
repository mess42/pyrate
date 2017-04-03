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
import numpy as np
import math

from globalconstants import standard_wavelength

class RayBundleNew(object):
    def __init__(self, x0, k0, Efield0, rayID = [], wave = standard_wavelength):
        """
        Class representing a bundle of rays.

        :param x0:  (2d numpy 3xN array of float) 
                    Initial position of the rays in global coordinates.  
        :param k0:  (2d numpy 3xN array of float) 
                    Wave vector of the rays in global coordinates.
        :param Efield0:  (2d numpy 3xN array of complex) 
                    Polarization state of the rays. 
        :param rayID: (1d numpy array of int) 
                    Set an ID number for each ray in the bundle; 
                    if empty -> generate arange 
        :param wave: (float) 
                    Wavelength of the radiation in millimeters. 
        """
        numray = np.shape(x0)[1]
        if rayID == [] or len(rayID) == 0:
            rayID = np.arange(numray)
        self.rayID = rayID            
            
        newshape = self.newshape(np.shape(x0))
            
        self.x = x0.reshape(newshape) 
        # First index counting index: x[0] == x0
        # shape(x): axis=0: counting axis
        # axis=1: vector components (xyz)
        # axis=2: ray number

        self.k = k0.reshape(newshape) 
        
        self.valid = np.ones((1, numray), dtype=bool)
        
        self.wave = wave
        if Efield0 == [] or len(Efield0) == 0:
            self.Efield = np.zeros(newshape)
            self.Efield[:, 1, :] = 1.
        else:
            self.Efield = Efield0.reshape(newshape)

    def newshape(self, shape2d):
        """
        Constructs 3d array shape (1, N, M) from 2d array shape (N, M).
        """
        return tuple([1] + list(shape2d))

    def append(self, xnew, knew, Enew, Validnew):
        """
        Appends one point with appropriate wave vector, electrical field and
        validity array. New validity status is cumulative. 
        
        :param xnew (2d numpy 3xN array of float)
        :param knew (2d numpy 3xN array of complex)
        :param Enew (2d numpy 3xN array of complex)
        :param Validnew (1d numpy array of bool)
        
        """
        newshape = self.newshape(np.shape(xnew))        
        newshapev = (1, np.shape(Validnew)[0])                
        
        txnew = np.reshape(xnew, newshape)
        tknew = np.reshape(knew, newshape)
        tEnew = np.reshape(Enew, newshape)
        tValidnew = np.reshape(self.valid[-1]*Validnew, newshapev)
        
        self.x      = np.vstack((self.x, txnew))
        self.k      = np.vstack((self.k, tknew))
        self.Efield = np.vstack((self.Efield, tEnew))
        self.valid  = np.vstack((self.valid, tValidnew))
        
    def clone(self):
        result = RayBundleNew(self.x[0], self.k[0], self.Efield[0], self.rayID, self.wave)
        
        result.x = np.copy(self.x)
        result.k = np.copy(self.k)
        result.Efield = np.copy(self.Efield)
        result.valid = np.copy(self.valid)

        return result        
        
        
    def getArcLength(self):
        """
        Calculates arc length for all rays.
        
        :return Arc length (1d numpy array of float)
        """
        ds = np.sqrt(np.sum((self.x[1:] - self.x[:-1])**2, axis=1)) # arc element lengths for every ray
        return np.sum(ds, axis=0)
        
    def returnKtoD(self):
        # S_j = Re((conj(E)_i E_i delta_{jl} - conj(E)_j E_l) k_l)
        absE2 = np.sum(np.conj(self.Efield)*self.Efield, axis=1)
        Ek = np.sum(self.Efield * self.k, axis=1)
        S = np.zeros_like( self.k, dtype=float)        
        for j in np.arange(3):
            S[:,j,:] = np.real(absE2*self.k[:,j,:] - Ek*np.conj(self.Efield)[:,j,:])
        absS = np.sqrt(np.sum(S**2, axis=1))
        d = np.zeros_like( self.k, dtype=float)        
        for j in np.arange(3):
            d[:,j,:] = S[:,j,:] / absS
        return d
        
    def draw2d(self, ax, color="blue", plane_normal = np.array([1, 0, 0]), up = np.array([0, 1, 0])):

        # normalizing plane_normal
        plane_normal = plane_normal/np.linalg.norm(plane_normal)
        up = up/np.linalg.norm(up)

        ez = np.cross(plane_normal, up)

        (num_points, num_dims, num_rays) = np.shape(self.x)

        # arrange num_ray copies of simple vectors in appropriate form
        plane_normal = np.column_stack((plane_normal for i in np.arange(num_rays)))
        ez = np.column_stack((ez for i in np.arange(num_rays)))
        up = np.column_stack((up for i in np.arange(num_rays)))
        
        ptlist = [self.x[i] for i in np.arange(num_points)]
        
        for (pt1, pt2) in zip(ptlist[1:], ptlist[:-1]):

            # perform in-plane projection
            pt1inplane = pt1 - np.sum(pt1*plane_normal, axis=0)*plane_normal
            pt2inplane = pt2 - np.sum(pt2*plane_normal, axis=0)*plane_normal
            
            # calculate y-components
            ypt1 = np.sum(pt1inplane * up, axis=0)
            ypt2 = np.sum(pt2inplane * up, axis=0)

            # calculate z-components
            zpt1 = np.sum(pt1inplane * ez, axis=0)
            zpt2 = np.sum(pt2inplane * ez, axis=0)
            
            y = np.vstack((ypt1, ypt2))
            z = np.vstack((zpt1, zpt2))
            ax.plot(z, y, color)            
            
        
class RayBundle(object):
    def __init__(self, o, d, mat, rayID, wave=standard_wavelength, pol=[]):
        """
        Class representing a bundle of rays.

        :param o:     Origin of the rays.  (2d numpy 3xN array of float)
        :param d:     Direction of the rays, normalized. (2d numpy 3xN array of float)
                      Direction of energy transport.
        :param rayID: Set an ID number for each ray in the bundle (1d numpy array of int)
                      (for example, the ray index at surface 0)
        :param wave:  Wavelength of the radiation in millimeters. (float)
        :param pol:   Polarization state of the rays. (2d numpy 2xN array of complex); not implemented yet

        """
        # Primary goal:
        # TODO: new properties: x0, k0; for GRIN media many xi, ki; last xN-1, kN-1        
        # TODO: remove t and calculate t by calcArcLength() which is a line integral t = int[x0, xN-1] ds

        # Secondary goal:
        # TODO: implement polarization
        # TODO: implement reflection / transmission coefficients
        #       coherent (Jones formalism) or incoherent (Mueller/Stokes formalism) ?
        #       the difference is, that Jones allows for description of phases of the coherent, fully polarized beam
        #       and Mueller allows for describing intensity transmission of partially polarized beams
        
        self.o = o
        self.d = d
        self.k = mat.returnDtoK(d, wave)
        self.rayID = rayID
        self.t = zeros(shape(o)[1])  # Geometrical path length to the ray final position.
        self.wave = wave
        self.pol = pol

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

    def draw2d(self, ax, color="blue"):
        # o and k in global coordinates
        nrays = shape(self.o)[1]
        for i in arange(nrays):
            y = array([self.o[1, i], self.o[1, i] + self.t[i] * self.d[1, i]])
            z = array([self.o[2, i], self.o[2, i] + self.t[i] * self.d[2, i]])
            ax.plot(z, y, color)


class RayPathNew(object):
    
    def __init__(self, initialraybundle):
        self.raybundles = [initialraybundle]
        
    def appendRayBundle(self, raybundle):
        self.raybundles.append(raybundle)
        
    def appendRayPath(self, raypath):
        self.raybundles += raypath.raybundles
        
    def draw2d(self, ax, color="blue", plane_normal = np.array([1, 0, 0]), up = np.array([0, 1, 0])):
        for r in self.raybundles:
            r.draw2d(ax, color=color, plane_normal=plane_normal, up=up)

        
        

class RayPath(object):
    def __init__(self, initialraybundle, opticalSystem):
        """
        Class representing the Path of a RayBundle through the whole optical system.

        :param initialraybundle: Raybundle at initial position in the optical system ( RayBundle object )
        :param opticalSystem:  optical system through which the rays are propagated ( OpticalSystem object )

        """
        self.raybundles = opticalSystem.trace(initialraybundle)

        #self.raybundles = [initialraybundle]
        #N = opticalSystem.getNumberOfSurfaces()
        #for i in arange(N-1)+1:
        #    self.traceToNextSurface(opticalSystem.surfaces[i-1], opticalSystem.surfaces[i])
           

    def traceToNextSurface(self, actualSurface, nextSurface):
        """
        Private routine that propagates a ray bundle to the next surface.
        Should call material.propagator from actualSurface.
        Please respect the privacy of this class and call it only from methods inside this class.
        intersection and normal are calculated in global coordinates.

        :param actualSurface: (Surface object)
        :param nextSurface: (Surface object)
        """

        intersection, t, normal, validIndices = \
                actualSurface.material.propagate(actualSurface, \
                                                nextSurface, \
                                                self.raybundles[-1])

        self.raybundles.append(nextSurface.material.refract(actualSurface.material, self.raybundles[-1], intersection, normal, validIndices))

    def draw2d(self, opticalsystem, ax, color="blue"):
        Nsurf = len(self.raybundles)
        for i in arange(Nsurf):
            self.raybundles[i].draw2d(ax, color=color)

if __name__ == "__main__":
    wavelength = standard_wavelength
    nray = 4
    x0      =       np.random.random((3,nray))
    k0      = 0.5 * np.random.random((3,nray))
    k0[2,:] = 1 - np.sqrt( k0[0,:]**2 + k0[1,:]**2 )
    k0      = k0 * 2 *pi / wavelength    
    
    E0 = np.random.random((3,nray)) # warning: E not orthogonal to k

    x1 = np.random.random((3,nray))
    
    validity = np.ones(nray, dtype=bool)

    r = RayBundleNew(x0, k0, E0, wave=wavelength)
    r.append(x1, k0, E0, validity)
    
    d = r.returnKtoD()
    print "d=",r.returnKtoD()
    print np.sum(d**2, axis=1)

    print(x0)
    print(x1)
    print(r.getArcLength())    
    
    
    r2 = RayBundleNew(x0 = [], k0 = [], Efield0 = [])    