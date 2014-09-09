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

import shape as surfShape # the name 'shape' already denotes the dimensions of a numpy array
import material
import pupil
import inspector

from numpy import *
from optimize import ClassWithOptimizableVariables

class Surface(ClassWithOptimizableVariables):
    """
    Represents a surface of an optical system.
    
    :param shap: Shape of the surface. Claculates the intersection with rays. ( Shape object or child )
    :param mater: Material of the volume behind the surface. Claculates the refraction. ( Material object or child )
    :param thickness: distance to next surface on the optical axis
    """
    def __init__(self, thickness = 0.0):
        self.listOfOptimizableVariables = []

        self.shap  = self.setShape("Conic")
        self.mater = self.setMaterial("ConstantIndexGlass")
        self.thickness = self.createOptimizableVariable("thickness", value = thickness, status=False)      

    def setThickness(self, thickness):
        self.thickness.val = thickness

    def getThickness(self):
        return self.thickness.val

    def setMaterial(self, materialType):
        """
        Sets the material object self.mater
        
        :param materialType: name of the material dispersion formula (str)

        :return self.mater: new Material object 
        """

       # conserve the most basic parameters of the shape
        try:        
            varsToRemove = self.mater.getAllOptimizableVariables()
            for v in varsToRemove:
                self.listOfOptimizableVariables.remove(v)
        except:
            pass

        names, classes = inspector.getListOfClasses(material, "<class \'material.", "<class \'material.Material\'>")

        self.mater = inspector.createObjectFromList(names, classes, materialType )

        # add optimizable variables of new shape
        self.listOfOptimizableVariables += self.mater.getAllOptimizableVariables()

        return self.mater

    def setMaterialCoefficients(self, coeff):
        """
        Sets the coefficients that determine the material behavior. 

        :param coeff: coefficients. Type and format depend on Material child class.
        """
        self.mater.setCoefficients(coeff)

    def setShape(self, shapeName):
        """
        Sets the shape object self.shap
        
        :param shapeName: name of the shape type (str)
        
        :return self.shap: new Shape object
        """

        # conserve the most basic parameters of the shape
        try:
            curv = self.shap.curvature.val 
            semidiam = self.shap.sdia.val
            
            varsToRemove = self.shap.getAllOptimizableVariables()
            for v in varsToRemove:
                self.listOfOptimizableVariables.remove(v)

        except:
            # self.shap does not exist yet
            curv = 0.0
            semidiam = 0.0 

        names, classes = inspector.getListOfClasses(surfShape, "<class \'shape.", "<class \'shape.Shape\'>")
 
        self.shap = inspector.createObjectFromList(names, classes, shapeName )
        self.shap.curvature.val = curv
        self.shap.sdia.val = semidiam

        # add optimizable variables of new shape
        self.listOfOptimizableVariables += self.shap.getAllOptimizableVariables()

        return self.shap

    def draw2d(self, ax, offset = [0,0], vertices=100, color="grey"):
        self.shap.draw2d(ax, offset, vertices, color)      

    def getABCDMatrix(self, nextSurface, ray):
        """
        Returns an ABCD matrix of the current surface.
        The matrix is set up in geometric convention for (y, dy/dz) vectors.

        The matrix contains:
        - paraxial refraction from vacuum through the front surface
        - paraxial translation through the material
        - paraxial refraction at the rear surface into vacuum

        :param nextSurface: next surface for rear surface curvature (Surface object)
        :param ray: ray bundle to obtain wavelength (RayBundle object)
        :return abcd: ABCD matrix (2d numpy 2x2 matrix of float)
        """
        curvature = self.shap.getCentralCurvature()
        nextCurvature = nextSurface.shap.getCentralCurvature()
        return self.mater.getABCDMatrix(curvature, self.thickness.val, nextCurvature, ray)
   
class OpticalSystem(ClassWithOptimizableVariables):
    """
    Represents an optical system, consisting of several surfaces and materials inbetween.
    """
    def __init__(self):
        self.surfaces = []
        self.insertSurface(0) # object
        self.insertSurface(1) # image
        self.surfaces[1].shap.sdia.val = 1E100
 
    def insertSurface(self,position):
        """
        Inserts a new surface into the optical system.

        :param position: number of the new surface (int). 
           Surface that is currently at this position 
           and all following surface indices are incremented.
        """
        self.surfaces.insert(position, Surface() )
        
    def removeSurface(self,position):
        """
        Removes a surface from the optical system.

        :param position: number of the surface to remove (int) 
        """
        self.surfaces.pop(position)
        
    def getNumberOfSurfaces(self):
        """
        Returns the number of surfaces, including object and image (int)
        """
        return len(self.surfaces)
        
    def setThickness(self, position, thickness):
        """
        Sets the on-axis thickness of a surface.

        :param position: number of the surface (int) 
        """
        self.surfaces[position].setThickness(thickness)

    def setMaterial(self, position, materialType):
        """
        Sets the material of a surface.

        :param position: number of the surface (int)
        :param materialType: name of the Material child class (str)
        """
        self.surfaces[position].setMaterial(materialType)

    def setMaterialCoefficients(self, position, coeff):
        """
        Sets the coefficients that determine the material behavior. 

        :param position: number of the surface (int)
        :param coeff: coefficients. Type and format depend on Material child class.
        """
        self.surfaces[position].setMaterialCoefficients(coeff)

    def setShape(self, position, shapeName):
        """
        Sets the shape of a surface.

        :param position: number of the surface (int)
        :param shapeName: name of the Shape child class (str)
        """
        self.surfaces[position].setShape(shapeName)

    def getABCDMatrix(self, ray, firstSurfacePosition = 0, lastSurfacePosition = -1):
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

        abcd = [[1,0],[0,1]]

        for i in arange( lastSurfacePosition - firstSurfacePosition + 1) + firstSurfacePosition:
            abcd = dot( self.surfaces[i].getABCDMatrix(self.surfaces[i+1], ray)  ,  abcd )

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
        abcdObjStop = self.getABCDMatrix(ray, 0 , stopPosition - 1) # object to stop

        zen  = abcdObjStop[0,1] / abcdObjStop[0,0] # entrance pupil position from object
        magen = 1.0 / abcdObjStop[0,0]     

        abcdStopIm = self.getABCDMatrix(ray, stopPosition, -1) # stop to image

        zex = - abcdStopIm[0,1] / abcdStopIm[1,1] # exit pupil position from image
        magex = abcdStopIm[0,0] - abcdStopIm[0,1] * abcdStopIm[1,0] / abcdStopIm[1,1]

        return zen, magen, zex, magex, abcdObjStop, abcdStopIm

    def getEffectiveFocalLength(self, ray):
        """
        Returns the effective (paraxial) focal length of the system.

        :param ray: Raybundle object
        :return f: focal length (float)
        """
        abcd = self.getABCDMatrix(ray)
        return -1.0 / abcd[1,0]

    def getParaxialMagnification(self, ray):
        """
        Returns the paraxial real space magnification of the system.
        Before calculation, the image is shifted into paraxial   finite conjugate plane.
 
        :param ray: Raybundle object
        :return pmag: real space paraxial magnification (float)
        """
        abcd = self.getABCDMatrix(ray)
        print abcd
        return abcd[0,0] - abcd[0,1] * abcd[1,0] / abcd[1,1]

    def draw2d(self, ax, offset = [0,0], vertices=100, color="grey"):
        N = self.getNumberOfSurfaces()
        offy = offset[0]
        offz = offset[1]
        for i in arange(N-1):
            self.surfaces[i].draw2d(ax, offset = [offy, offz])
            offz += self.surfaces[i].getThickness()
 

    def createOptimizableVariable(self, name, value = 0.0, status=False):
        """
        This class is not able to create own variables. 
        It only forwards variables from its surfaces.
        """
        raise NotImplementedError()
     
    def getAllOptimizableVariables(self):
        varsToReturn = []
        for sur in self.surfaces:
            varsToReturn += sur.getAllOptimizableVariables()
        return varsToReturn
        



