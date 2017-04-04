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
from material import ConstantIndexGlass

import aperture
import pupil
from coordinates import LocalCoordinates

#import inspector
import numpy as np
from optimize import ClassWithOptimizableVariables
from optimize import OptimizableVariable

from ray import RayPathNew, RayBundleNew

from globalconstants import canonical_ex, canonical_ey, canonical_ez

import uuid


class CoordinateTreeBase(ClassWithOptimizableVariables):
    """
    Optical element base class for optical system and optical element, whereas
    an optical system consists of many optical elements.
    Implements functionality for local coordinate system tree and connection
    checks.
    
    :param rootcoordinatesystem (LocalCoordinates object)
    :param label (string)
    :param *kwargs (key word arguments)
    """
    def __init__(self, rootcoordinatesystem, label="", **kwargs):
        self.label = label
        self.rootcoordinatesystem = rootcoordinatesystem
        
    def checkForRootConnection(self, lc):
        """
        Checks whether given local coordinate system is child of rootcoordinatesystem.
        
        :param lc (LocalCoordinates object)
        
        :return bool
        """
        allconnectedchildren = self.rootcoordinatesystem.returnConnectedChildren()        
        return (lc in allconnectedchildren)
            
    def addLocalCoordinateSystem(self, lc, refname):
        """
        Adds local coordinate system as child to given reference.
        
        :param lc (LocalCoordinates object)
        :param refname (string)
        
        :return lc        
        """
        allnames = self.rootcoordinatesystem.returnConnectedNames()
       
        if lc.name in allnames:
            lc.name = str(uuid.uuid4())
            
        if refname not in allnames:
            refname = self.rootcoordinatesystem.name
        
        self.rootcoordinatesystem.addChildToReference(refname, lc)
            
        return lc
    
    
    
        


class OpticalElement(CoordinateTreeBase):
    """
    Represents an optical element (volume with surface boundary and inner
    surfaces representing material boundaries)
    
    :param lc (Local Coordinates of optical element)
    :param label (string), if empty -> uuid
    """
    def __init__(self, lc, label="", **kwargs):
        super(OpticalElement, self).__init__(lc, label=label)
        self.__surfaces = {} # Append surfaces objects
        self.__materials = {} # Append materials objects
        self.__surf_mat_connection = {} # dict["surfname"] = ("mat_minus_normal", "mat_plus_normal")
    
    
    def addSurface(self, key, surface_object, (minusNmat_key, plusNmat_key), label=""):
        """
        Adds surface class object to the optical element.
        
        :param key (string ... dict key)
        :param surface_object (Surface class object)
        :param (minusNmat_key, plusNmat_key) (tuple of strings ... keys to material dict)
        :param label (string, optional), label of surface
        """
        if self.checkForRootConnection(surface_object.rootcoordinatesystem):
            self.__surfaces[key] = surface_object
        else:
            raise Exception("surface coordinate system should be connected to OpticalElement root coordinate system")
        self.__surfaces[key].label = label
        self.__surf_mat_connection[key] = (minusNmat_key, plusNmat_key)

    def getSurfaces(self):
        return self.__surfaces
        
    surfaces = property(fget=getSurfaces)
        

    def addMaterial(self, key, material_object, comment=""):
        """
        Adds material class object to the optical element.

        :param key (string ... dict key)
        :param material_object (Material class object)
        :param comment (string, optional), comment for the material
        """
        if self.checkForRootConnection(material_object.lc):
            self.__materials[key] = material_object
        else:
            raise Exception("material coordinate system should be connected to OpticalElement root coordinate system")            
        self.__materials[key].comment = comment

    def findoutWhichMaterial(self, mat1, mat2, current_mat):
        """
        Dirty method to determine material after refraction. 
        (Reference comparison.)
        
        :param mat1 (Material object)
        :param mat2 (Material object)
        :param current_mat (Material object)
        
        :return (Material object)
        """        
        
        if id(mat1) == id(current_mat):
            returnmat = mat2
        else:
            returnmat = mat1
            
        return returnmat

     

    def seqtrace(self, raybundle, sequence, background_medium):
        
        current_material = background_medium    
    
        rpath = RayPathNew(raybundle)    
    
        for (surfkey, refract_flag, ordinary_flag) in sequence:
            
            current_bundle = rpath.raybundles[-1]
            current_surface = self.__surfaces[surfkey]
            
            #print(rpath.raybundles[-1].x[-1, :, 0].reshape((3, 1)))
            #print(current_surface.shape.getNormalDerivative(rpath.raybundles[-1].x[-1, :, 0].reshape((3, 1))))
            
            current_material.propagate(current_bundle, current_surface)
            
            (mnmat, pnmat) = self.__surf_mat_connection[surfkey]
            mnmat = self.__materials.get(mnmat, background_medium)
            pnmat = self.__materials.get(pnmat, background_medium)

            if refract_flag:
                current_material = self.findoutWhichMaterial(mnmat, pnmat, current_material)
                rpath.appendRayBundle(current_material.refractNew(current_bundle, current_surface))
            else:
                rpath.appendRayBundle(current_material.reflectNew(current_bundle, current_surface))
                


        return rpath



class SurfaceNew(CoordinateTreeBase):
    """
    Represents a surface of an optical system.

    :param shape: Shape of the surface. Calculates the intersection with rays. ( Shape object or child )
    :param material: Material of the volume behind the surface. Calculates the refraction. ( Material object or child )
    :param thickness: distance to next surface on the optical axis
    """
    def __init__(self, rootlc, shape=None, apert=None, label=""):
        super(SurfaceNew, self).__init__(rootlc, label)

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
        
        self.shape.intersectNew(raybundle)
        
        if remove_rays_outside_aperture:
            globalintersection = raybundle.x[-1]
            local_ap_intersection = self.aperture.lc.returnGlobalToLocalPoints(globalintersection)
        
            valid = self.aperture.arePointsInAperture(local_ap_intersection[0], local_ap_intersection[1])
        
            raybundle.valid[-1] = raybundle.valid[-1]*valid


    def draw2d(self, ax, vertices=100, inyzplane = True, color="grey", plane_normal = canonical_ex, up = canonical_ey):
        """
        :param ax (Axis object)
        :param vertices (int), vertices in xy for aperture sampling
        :param inyzplane (bool), cuts globalpts in yz plane before projection on plane_normal
        :param color (string), "red", "blue", "grey", "green", ...
        :param plane_normal (1D numpy array of float), new x projection axis
        :param up (1D numpy array), invariant y axis, z = x x y
        
        """


        sizelimit = 1000.0
        failsafevalue = 10.0        
        if self.aperture == None:
            effsemidia = failsafevalue
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
        

class OpticalSystemNew(CoordinateTreeBase):
    """
    Represents an optical system, consisting of several surfaces and materials inbetween.
    """
    def __init__(self, rootlc = None, matbackground = None, name = ""):
        # TODO: rename variable name to label
        """
        Creates an optical system object. Initially, it contains 2 plane surfaces (object and image).


        :param objectLC: local coordinate system of object (class LocalCoordinates).


        """
        if rootlc is None:        
            rootlc = LocalCoordinates(name="global")
        self.rootcoordinatesystem = rootlc

        super(OpticalSystemNew, self).__init__(self.rootcoordinatesystem, label = name)
        
        if matbackground is None:
            matbackground = ConstantIndexGlass(self.rootcoordinatesystem, 1.0)

        self.material_background = matbackground # Background material        
        self.elements = {}

    def seqtrace(self, initialbundle, elementsequence): # [("elem1", [1, 3, 4]), ("elem2", [1,4,4]), ("elem1", [4, 3, 1])]
        rpath = RayPathNew(initialbundle)
        for (elem, subseq) in elementsequence:
            rpath.appendRayPath(self.elements[elem].seqtrace(rpath.raybundles[-1], subseq, self.material_background)) 
        return rpath
            

    def addElement(self, key, element):
        """
        Adds a new element (containing several surfaces) into the optical system.

        :param key (string)        
        :param element (optical element class)
        """
        if self.checkForRootConnection(element.rootcoordinatesystem):
            self.elements[key] = element
        else:
            raise Exception("OpticalElement root should be connected to root of OpticalSystem")

    def removeElement(self, key):
        """
        Removes an optical element from the optical system.

        :param key (string)
        """
        # TODO: update of local coordinate references missing
        if key in self.elements:        
            self.elements.pop(key)


    def getABCDMatrix(self, ray, firstSurfacePosition=0, lastSurfacePosition=-1):
        """
        Returns an ABCD matrix of the optical system.
        The matrix is set up in geometric convention for (y, dy/dz) vectors.

        The matrix contains:
        - paraxial refraction from vacuum through the first surface
        - paraxial propagation through the system
        - paraxial refraction after the last surface into vacuum

        :param firstSurfacePosition: Position of the first surface to consider (int).
          Preset is 0 (object position).
        :param lastSurfacePosition: Position of the last surface to consider (int).
          Preset is -1 (image position)
        :param ray: Ray bundle object.
        :return abcd: ABCD matrix (2d numpy 2x2 matrix of float)
        """

        if lastSurfacePosition < 0:
            lastSurfacePosition = self.getNumberOfSurfaces() - lastSurfacePosition - 3

        abcd = [[1, 0], [0, 1]]

        for i in np.arange(lastSurfacePosition - firstSurfacePosition + 1) + firstSurfacePosition:
            abcd = np.dot(self.surfaces[i].getABCDMatrix(self.surfaces[i+1], ray), abcd)

        return abcd



    def getParaxialPupil(self, stopPosition, ray):
        """
        Returns the paraxially calculated pupil positions.

        :param stopPosition: surface number of the surface that is defined as stop (int)
        :param ray: Raybundle object

        :return zen: entrance pupil position from object (float)
        :return magen: entrance pupil magnificaction; entrance pupil diameter per stop diameter (float)
        :return zex: exit pupil position from image (float)
        :return magex: exit pupil magnificaction; exit pupil diameter per stop diameter (float)
        """
        abcdObjStop = self.getABCDMatrix(ray, 0, stopPosition - 1)  # object to stop

        zen = abcdObjStop[0, 1] / abcdObjStop[0, 0]  # entrance pupil position from object
        magen = 1.0 / abcdObjStop[0, 0]

        abcdStopIm = self.getABCDMatrix(ray, stopPosition, -1)  # stop to image

        zex = - abcdStopIm[0, 1] / abcdStopIm[1, 1]  # exit pupil position from image
        magex = abcdStopIm[0, 0] - abcdStopIm[0, 1] * abcdStopIm[1, 0] / abcdStopIm[1, 1]

        return zen, magen, zex, magex, abcdObjStop, abcdStopIm

    def getEffectiveFocalLength(self, ray):
        """
        Returns the effective (paraxial) focal length of the system.

        :param ray: Raybundle object
        :return f: focal length (float)
        """
        abcd = self.getABCDMatrix(ray)
        return -1.0 / abcd[1, 0]

    def getParaxialMagnification(self, ray):
        """
        Returns the paraxial real space magnification of the system.
        Before calculation, the image is shifted into paraxial   finite conjugate plane.

        :param ray: Raybundle object
        :return pmag: real space paraxial magnification (float)
        """
        abcd = self.getABCDMatrix(ray)
        print abcd
        return abcd[0, 0] - abcd[0, 1] * abcd[1, 0] / abcd[1, 1]


    def draw2d(self, ax, vertices=100, color="grey"):
        for (num, s) in enumerate(self.surfaces):
            s.draw2d(ax, vertices=vertices, color=color)
            


if __name__ == "__main__":

    # AC254-100-Ad 	25.4 	100.1 	97.1 	info 	62.8 	-45.7 	-128.2 	4.0 	2.5 	4.7 	N-BK7/SF5    
    
    os = OpticalSystemNew()
    
    lc1 = os.addLocalCoordinateSystem(LocalCoordinates(decz=10.0), refname=os.rootcoordinatesystem.name)
    lc2 = os.addLocalCoordinateSystem(LocalCoordinates(decz=20.0), refname=lc1.name)
    lc3 = os.addLocalCoordinateSystem(LocalCoordinates(decz=30.0), refname=lc2.name)
    lc4 = os.addLocalCoordinateSystem(LocalCoordinates(decz=40.0), refname=lc3.name)
    
    lc5 = os.addLocalCoordinateSystem(LocalCoordinates(name="COM", decx=10.0, decy=5.0, decz=10.), refname=lc1.name)
    
    print(os.rootcoordinatesystem.pprint())
        
