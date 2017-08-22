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

import surfShape
import aperture
from localcoordinatestreebase import LocalCoordinatesTreeBase
import numpy as np
from globalconstants import canonical_ex, canonical_ey



class Surface(LocalCoordinatesTreeBase):
    """
    Represents a surface of an optical system.

    :param shape: Shape of the surface. Calculates the intersection with rays. ( Shape object or child )
    :param material: Material of the volume behind the surface. Calculates the refraction. ( Material object or child )
    :param thickness: distance to next surface on the optical axis
    """
    def __init__(self, rootlc, shape=None, apert=None, **kwargs):
        super(Surface, self).__init__(rootlc, **kwargs)

        if shape is None:        
            shape = surfShape.Conic(rootlc)
        if apert is None:
            apert = aperture.BaseAperture(rootlc)
            
        self.setShape(shape)
        self.setAperture(apert)


    def setAperture(self, apert):
        """
        Sets the shape object self.shap

        :param shape: the new Shape object

        :return self.shape: new Shape object
        """
        if self.checkForRootConnection(apert.lc):
            self.__aperture = apert
        else:
            raise Exception("Aperture coordinate system should be connected to surface coordinate system")            
       
    def getAperture(self):
        return self.__aperture
        
    aperture = property(getAperture, setAperture)


    def setShape(self, shape):
        """
        Sets the shape object self.shap

        :param shape: the new Shape object

        :return self.shape: new Shape object
        """
        if self.checkForRootConnection(shape.lc):
            self.__shape = shape
        else:
            raise Exception("Shape coordinate system should be connected to surface coordinate system")            
        
    def getShape(self):
        return self.__shape
        
    shape = property(getShape, setShape)
    

    def intersect(self, raybundle, remove_rays_outside_aperture=True):
        """
        Calculates intersection from raybundle. Knows shape and aperture and
        can remove rays due to aperture.
        
        :param raybundle (RayBundle object), gets changed!
        """
        
        self.shape.intersect(raybundle)
        
        if remove_rays_outside_aperture:
            globalintersection = raybundle.x[-1]
            local_ap_intersection = self.aperture.lc.returnGlobalToLocalPoints(globalintersection)
        
            valid = self.aperture.arePointsInAperture(local_ap_intersection[0], local_ap_intersection[1])
        
            raybundle.valid[-1] = raybundle.valid[-1]*valid


    def draw2d(self, ax, vertices=100, inyzplane=True, color="grey", plane_normal=canonical_ex, up=canonical_ey):
        """
        :param ax (Axis object)
        :param vertices (int), vertices in xy for aperture sampling
        :param inyzplane (bool), cuts globalpts in yz plane before projection on plane_normal
        :param color (string), "red", "blue", "grey", "green", ...
        :param plane_normal (1D numpy array of float), new x projection axis
        :param up (1D numpy array), invariant y axis, z = x x y
        
        """


        sizelimit = 1000.0
        failsafevalue = 11.0        
        if self.aperture == None:
            effsemidia = failsafevalue
            # TODO: choose max ray height of all bundles instead 
            # ( cosmetic but absolutely necessary for beauty )
        else:
            if self.aperture.getTypicalDimension() <= sizelimit:
                # TODO: maybe introduce aperture types Object and Image to distuingish from very large normal apertures
                effsemidia = self.aperture.getTypicalDimension() #self.sdia.val if self.sdia.val < 10.0 else 10.0
            else:
                effsemidia = failsafevalue
        
        xl = effsemidia * np.linspace(-1, 1, num=vertices)
        yl = effsemidia * np.linspace(-1, 1, num=vertices)
        
        X, Y = np.meshgrid(xl, yl)
        x = X.flatten()
        y = Y.flatten()
        
        isinap = np.array(self.aperture.arePointsInAperture(x, y))
        xinap = x[isinap]        
        yinap = y[isinap]
        zinap = np.zeros_like(xinap)
        
        localpts_aperture = np.row_stack((xinap, yinap, zinap))
        localpts_shape = self.shape.lc.returnOtherToActualPoints(localpts_aperture, self.aperture.lc)
        
        xinap_shape = localpts_shape[0, :]
        yinap_shape = localpts_shape[1, :]        
        zinap_shape = self.shape.getSag(xinap_shape, yinap_shape)
        
        localpts_shape = np.row_stack((xinap_shape, yinap_shape, zinap_shape))
        localpts_surf = self.rootcoordinatesystem.returnOtherToActualPoints(localpts_shape, self.shape.lc)        
        
        # ebenenprojektion hier!        
        
        globalpts = self.rootcoordinatesystem.returnLocalToGlobalPoints(localpts_surf)

        # doubled code begin (also in RayBundleNew.draw2d)
        plane_normal = plane_normal/np.linalg.norm(plane_normal)
        up = up/np.linalg.norm(up)

        ez = np.cross(plane_normal, up)

        (num_dims, num_rays) = np.shape(globalpts)

        # arrange num_ray copies of simple vectors in appropriate form
        plane_normal = np.column_stack((plane_normal for i in np.arange(num_rays)))
        ez = np.column_stack((ez for i in np.arange(num_rays)))
        up = np.column_stack((up for i in np.arange(num_rays)))
        # doubled code (also in RayBundleNew.draw2d)

        # doubled code (see ray.py)

        # show surfaces with points in yz plane before pn-plane projection
        # if false show full surface
        if inyzplane:
            inYZplane = np.abs(globalpts[0]) < 2*effsemidia/vertices
            globalpts = globalpts[:, inYZplane]
            plane_normal = plane_normal[:, inYZplane]
            up = up[:, inYZplane]
            ez = ez[:, inYZplane]

        globalptsinplane = globalpts - np.sum(globalpts*plane_normal, axis=0)*plane_normal

       
        # calculate y-components
        ypt = np.sum(globalptsinplane * up, axis=0)
        # calculate z-components
        zpt = np.sum(globalptsinplane * ez, axis=0)

        # doubled code (see ray.py)

        
        #ax.plot(zinap+offset[1], yinap+offset[0], color)
        ax.plot(zpt, ypt, color)
        
        
        #self.shape.draw2d(ax, offset, vertices, color, self.aperture)

    def getCentralCurvature(self, ray):
        curvature = self.shape.getCentralCurvature()
        # TODO: curvature at ray position
        
        return curvature
        

