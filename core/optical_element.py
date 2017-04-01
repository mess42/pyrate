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

import uuid

class OpticalElement(ClassWithOptimizableVariables):
    """
    Represents an optical element (volume with surface boundary and inner
    surfaces representing material boundaries)
    
    :param lc (Local Coordinates of optical element)
    :param label (string), if empty -> uuid
    """
    def __init__(self, lc, label="", **kwargs):
        self.label = label
        self.__surfaces = {} # Append surfaces objects
        self.__materials = {} # Append materials objects
        self.__surf_mat_connection = {} # dict["surfname"] = ("mat_minus_normal", "mat_plus_normal")
        
        self.lc = lc
        
    def addSurface(self, key, surface_object, (minusNmat_key, plusNmat_key), label=""):
        """
        Adds surface class object to the optical element.
        
        :param key (string ... dict key)
        :param surface_object (Surface class object)
        :param (minusNmat_key, plusNmat_key) (tuple of strings ... keys to material dict)
        :param label (string, optional), label of surface
        """
        self.__surfaces[key] = surface_object
        self.__surfaces[key].label = label
        self.__surf_mat_connection[key] = (minusNmat_key, plusNmat_key)
        

    def addMaterial(self, key, material_object, comment=""):
        """
        Adds material class object to the optical element.

        :param key (string ... dict key)
        :param material_object (Material class object)
        :param comment (string, optional), comment for the material
        """
        self.__materials[key] = material_object
        self.__materials[key].comment = comment

    def seqtrace(self, raybundle, sequence, background_medium):
        # TODO: hier weitermachen
        # sequence = ["surf1", "surf2", "surf3"], keys
    
        current_material = background_medium    
    
        for surfkey in sequence:
            current_material.propagate(self.__surfaces[surfkey])            
            print(surfkey)
            print(self.__surf_mat_connection[surfkey])
            (mnmat, pnmat) = self.__surf_mat_connection[surfkey]
            mnmat = self.__materials.get(mnmat, background_medium)
            pnmat = self.__materials.get(pnmat, background_medium)
            
            print(mnmat.n.evaluate(), pnmat.n.evaluate())
    
        return RayPathNew(raybundle)



class SurfaceNew(ClassWithOptimizableVariables):
    """
    Represents a surface of an optical system.

    :param shape: Shape of the surface. Calculates the intersection with rays. ( Shape object or child )
    :param material: Material of the volume behind the surface. Calculates the refraction. ( Material object or child )
    :param thickness: distance to next surface on the optical axis
    """
    def __init__(self, lc, shape=surfShape.Conic(), aperture=aperture.BaseAperture(), **kwargs):
        super(SurfaceNew, self).__init__()

        self.shape = shape
        self.aperture = aperture
        self.lc = lc # reference to local coordinate system tree

    

    def setShape(self, shape):
        """
        Sets the shape object self.shap

        :param shape: the new Shape object

        :return self.shape: new Shape object
        """

        # TODO: conserve the most basic parameters of the shape

        self.shape = shape
        return self.shape

    def draw2d(self, ax, offset=(0, 0), vertices=100, color="grey"):
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
        
        
        zinap = self.shape.getSag(xinap, yinap)
        
        localpts = np.row_stack((xinap, yinap, zinap))
        globalpts = self.lc.returnLocalToGlobalPoints(localpts)

        inYZplane = np.abs(xinap) < 2*effsemidia/vertices

        globalpts = globalpts[:, inYZplane]

        
        #ax.plot(zinap+offset[1], yinap+offset[0], color)
        ax.plot(globalpts[2], globalpts[1], color)
        
        
        #self.shape.draw2d(ax, offset, vertices, color, self.aperture)

    def getCentralCurvature(self, ray):
        curvature = self.shape.getCentralCurvature()
        # TODO: curvature at ray position
        
        return curvature
        

class OpticalSystemNew(ClassWithOptimizableVariables):
    """
    Represents an optical system, consisting of several surfaces and materials inbetween.
    """
    def __init__(self, matbackground = ConstantIndexGlass(1.0), name = "", objectLC = LocalCoordinates(name="object")):
        """
        Creates an optical system object. Initially, it contains 2 plane surfaces (object and image).


        :param objectLC: local coordinate system of object (class LocalCoordinates).


        """
        super(OpticalSystemNew, self).__init__(name = name)
        
        self.globalcoordinatesystem = LocalCoordinates(name="global")
        self.lcfocus = "global"
        
        self.objectlc = self.addLocalCoordinateSystem(objectLC)
        self.lcfocus = self.objectlc.name
        

        self.material_background = matbackground # Background material        
        self.elements = {}
        self.addElement("object", OpticalElement(self.objectlc))  # object
        # in standard initialization the surface use the BaseAperture which is not limited

    def seqtrace(self, initialbundle, elementsequence): # [("elem1", [1, 3, 4]), ("elem2", [1,4,4]), ("elem1", [4, 3, 1])]
        rpath = RayPathNew(initialbundle)
        for (elem, subseq) in elementsequence:
            rpath.appendRayPath(self.elements[elem].seqtrace(rpath.raybundles[-1], subseq, self.material_background)) 
        return rpath
            

    def addLocalCoordinateSystem(self, tmplc, refname=""):
        allnames = self.globalcoordinatesystem.returnConnectedNames()
       
        if refname == "":
            refname = self.lcfocus
        if tmplc.name in allnames:
            # TODO: throw exception
            tmplc.name = ""
            
        if refname not in allnames:
            refname = self.globalcoordinates.name
        
        self.globalcoordinatesystem.addChildToReference(refname, tmplc)
            
        self.lcfocus = tmplc.name
        
        return tmplc
            
    def addElement(self, key, element):
        """
        Adds a new element (containing several surfaces) into the optical system.

        :param key (string)        
        :param element (optical element class)
        """

        self.elements[key] = element

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
    
    lc1 = os.addLocalCoordinateSystem(LocalCoordinates(decz=10.0))
    os.addLocalCoordinateSystem(LocalCoordinates(decz=20.0))
    os.addLocalCoordinateSystem(LocalCoordinates(decz=30.0))
    os.addLocalCoordinateSystem(LocalCoordinates(decz=40.0))
    
    os.addLocalCoordinateSystem(LocalCoordinates(name="COM", decx=10.0, decy=5.0, decz=10.), refname=lc1.name)
    
    print(os.globalcoordinatesystem.pprint())
        
