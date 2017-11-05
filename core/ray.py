#!/usr/bin/env/python
"""
Pyrate - Optical raytracing based on Python

Copyright (C) 2014 Moritz Esslinger moritz.esslinger@web.de
               and Johannes Hartung j.hartung@gmx.net
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

from globalconstants import standard_wavelength, canonical_ex, canonical_ey

class RayBundle(object):
    def __init__(self, x0, k0, Efield0, rayID = [], wave = standard_wavelength, splitted=False):
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
        self.splitted = splitted
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
        if Efield0 is None or len(Efield0) == 0:
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
        result = RayBundle(self.x[0], self.k[0], self.Efield[0], self.rayID, self.wave)
        
        result.x = np.copy(self.x)
        result.k = np.copy(self.k)
        result.Efield = np.copy(self.Efield)
        result.valid = np.copy(self.valid)

        return result        
        
        
    def returnLocalComponents(self, lc, num):
        xloc = lc.returnGlobalToLocalPoints(self.x[num])
        kloc = lc.returnGlobalToLocalDirections(self.k[num])
        Eloc = lc.returnGlobalToLocalDirections(self.Efield[num])
        
        return (xloc, kloc, Eloc)

    def returnLocalD(self, lc, num):
        dloc = lc.returnGlobalToLocalDirections(self.returnKtoD()[num])
        return dloc
        
    def appendLocalComponents(self, lc, xloc, kloc, Eloc, valid):
        xglob = lc.returnLocalToGlobalPoints(xloc)
        kglob = lc.returnLocalToGlobalDirections(kloc)
        Eglob = lc.returnLocalToGlobalDirections(Eloc)
        
        self.append(xglob, kglob, Eglob, valid)
        
        
    def returnKtoD(self):
        # S_j = Re((conj(E)_i E_i delta_{jl} - conj(E)_j E_l) k_l)

        # FIXME: improve this by using np.einsum?
        
#        absE2 = np.sum(np.conj(self.Efield)*self.Efield, axis=1)
#        Ek = np.sum(self.Efield * self.k, axis=1)
#        S = np.zeros_like( self.k, dtype=float)
#        for j in np.arange(3):
#            S[:,j,:] = np.real(absE2*self.k[:,j,:] - Ek*np.conj(self.Efield)[:,j,:])
#        absS = np.sqrt(np.sum(S**2, axis=1))
#        d = np.zeros_like( self.k, dtype=float)        
#        for j in np.arange(3):
#            d[:,j,:] = S[:,j,:] / absS
#        return d

        absE2 = np.sum(np.conj(self.Efield)*self.Efield, axis=1)
        Ek = np.sum(self.Efield*self.k, axis=1)
        S = np.real(absE2*self.k - Ek*np.conj(self.Efield))
        return S / np.sqrt(np.sum(S**2, axis=1))

        
        
    def getLocalSurfaceNormal(self, surface, material, xglob):
        xlocshape = surface.shape.lc.returnGlobalToLocalPoints(xglob)
        nlocshape = surface.shape.getNormal(xlocshape[0], xlocshape[1])
        nlocmat = material.lc.returnOtherToActualDirections(nlocshape, surface.shape.lc)
        return nlocmat

        
    def draw2d(self, ax, color="blue", plane_normal=canonical_ex, up=canonical_ey):

        # normalizing plane_normal, up direction
        plane_normal = plane_normal/np.linalg.norm(plane_normal)
        up = up/np.linalg.norm(up)

        ez = np.cross(plane_normal, up)

        (num_points, num_dims, num_rays) = np.shape(self.x)
        
        if num_rays == 0:
            return

        # arrange num_ray copies of simple vectors in appropriate form
        plane_normal = np.column_stack((plane_normal for i in np.arange(num_rays)))
        ez = np.column_stack((ez for i in np.arange(num_rays)))
        up = np.column_stack((up for i in np.arange(num_rays)))
        
        ptlist = [self.x[i] for i in np.arange(num_points)]
        validity = [self.valid[i] for i in np.arange(num_points)]        
        
        for (pt1, pt2, todraw) in zip(ptlist[1:], ptlist[:-1], validity[1:]):

            # perform in-plane projection
            pt1inplane = pt1 - np.sum(pt1*plane_normal, axis=0)*plane_normal
            pt2inplane = pt2 - np.sum(pt2*plane_normal, axis=0)*plane_normal
            
            # calculate y-components
            ypt1 = np.sum(pt1inplane * up, axis=0)
            ypt2 = np.sum(pt2inplane * up, axis=0)

            # calculate z-components
            zpt1 = np.sum(pt1inplane * ez, axis=0)
            zpt2 = np.sum(pt2inplane * ez, axis=0)
            
            y = np.vstack((ypt1, ypt2))[:, todraw]
            z = np.vstack((zpt1, zpt2))[:, todraw]
            ax.plot(z, y, color)            
            
        


class RayPath(object):
    
    def __init__(self, initialraybundle=None):
        if initialraybundle is None:
            self.raybundles = []
        else:
            self.raybundles = [initialraybundle]
        
    def appendRayBundle(self, raybundle):
        self.raybundles.append(raybundle)
        
    def appendRayPath(self, raypath):
        self.raybundles += raypath.raybundles
        
    def draw2d(self, ax, color="blue", plane_normal=canonical_ex, up=canonical_ey):
        for r in self.raybundles:
            r.draw2d(ax, color=color, plane_normal=plane_normal, up=up)

    def containsSplitted(self):
        return any([r.splitted for r in self.raybundles])
        


if __name__ == "__main__":
    wavelength = standard_wavelength
    nray = 4
    x0      =       np.random.random((3,nray))
    k0      = 0.5 * np.random.random((3,nray))
    k0[2,:] = 1 - np.sqrt( k0[0,:]**2 + k0[1,:]**2 )
    k0      = k0 * 2 *math.pi / wavelength    
    
    E0 = np.random.random((3,nray)) # warning: E not orthogonal to k

    x1 = np.random.random((3,nray))
    
    validity = np.ones(nray, dtype=bool)

    r = RayBundle(x0, k0, E0, wave=wavelength)
    r.append(x1, k0, E0, validity)
    
    d = r.returnKtoD()
    print "d=",r.returnKtoD()
    print np.sum(d**2, axis=1)

    print(x0)
    print(x1)
    print(r.getArcLength())    
    
    
