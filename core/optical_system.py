#!/usr/bin/python3
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

import numpy as np
from optimize import ClassWithOptimizableVariables
from optimize import OptimizableVariable


class Surface(ClassWithOptimizableVariables):
    """
    Represents a surface of an optical system.

    :param shape: Shape of the surface. Calculates the intersection with rays. ( Shape object or child )
    :param material: Material of the volume behind the surface. Calculates the refraction. ( Material object or child )
    :param thickness: distance to next surface on the optical axis
    """
    def __init__(self, shape=surfShape.Conic(), thickness=0.0, material=ConstantIndexGlass(), aperture=aperture.BaseAperture()):
        super(Surface, self).__init__()

        self.shape = shape
        self.material = material
        self.aperture = aperture

        #self.thickness = self.createOptimizableVariable("thickness", value=thickness, status=False)
        #self.copyOptimizableVariables(shape)
        #self.copyOptimizableVariables(material)

        #self.thickness =
        self.addVariable("thickness", OptimizableVariable(False, "Variable", value=thickness))
        # TODO: new style code

    def setThickness(self, thickness):
        self.dict_variables["thickness"].setvalue(thickness)

    def getThickness(self):
        return self.dict_variables["thickness"].evaluate()

    def setMaterial(self, materialType):
        """
        Sets the material object self.mater

        :param materialType: name of the material dispersion formula (str)

        :return self.material: new Material object
        """

        # conserve the most basic parameters of the shape
        # TODO: should not be necessary anymore the old material will be overwritten and
        # the new materials dict will be appended to the variables list as necessary
        """
        try:
            varsToRemove = self.material.getAllOptimizableVariables()
            for v in varsToRemove:
                self.listOfOptimizableVariables.remove(v)
        except:
            pass

        self.material = materialType()

        print( "orig listofoptvars: ", [i.name for i in self.listOfOptimizableVariables] )
        print( "material listofoptvars: ", [i.name for i in self.material.getAllOptimizableVariables()] )

        # add optimizable variables of new shape
        self.listOfOptimizableVariables += self.material.getAllOptimizableVariables()

        print( "new listofoptvars: ", [i.name for i in self.listOfOptimizableVariables] )
        """

        return self.material

    def setMaterialCoefficients(self, coeff):
        """
        Sets the coefficients that determine the material behavior.

        :param coeff: coefficients. Type and format depend on Material child class.
        """
        self.material.setCoefficients(coeff)

    def setShape(self, shape):
        """
        Sets the shape object self.shap

        :param shape: the new Shape object

        :return self.shape: new Shape object
        """

        # conserve the most basic parameters of the shape
        """
        try:
            curv = self.shape.curvature.val
            #semidiam = self.shape.sdia.val

            varsToRemove = self.shape.getAllOptimizableVariables()
            for v in varsToRemove:
                self.listOfOptimizableVariables.remove(v)

        except:
            # self.shape does not exist yet
            curv = 0.0
            #semidiam = 0.0

        # names, classes = inspector.getListOfClasses(surfShape, "<class \'shape.", "<class \'shape.Shape\'>")

        # self.shape = inspector.createObjectFromList(names, classes, shapeName)
        """
        self.shape = shape
        # OLD
        #self.shape.curvature.val = curv


        #self.shape.sdia.val = semidiam


        # add optimizable variables of new shape

        # OLD
        #self.listOfOptimizableVariables += self.shape.getAllOptimizableVariables()

        return self.shape

    def draw2d(self, ax, offset=(0, 0), vertices=100, color="grey"):
        self.shape.draw2d(ax, offset, vertices, color, self.aperture)

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
        curvature = self.shape.getCentralCurvature()
        # TODO: improvement, call shape.getHessian() to obtain curvature components at intersection point
        nextCurvature = nextSurface.shape.getCentralCurvature()
        return self.material.getABCDMatrix(curvature, self.getThickness(), nextCurvature, ray)
        # TODO:


class OpticalSystem(ClassWithOptimizableVariables):
    """
    Represents an optical system, consisting of several surfaces and materials inbetween.
    """
    def __init__(self, objectDistance = 0.0, primaryWavelength = 550e-6):
        """
        Creates an optical system object. Initially, it contains 2 plane surfaces (object and image).

        :param objectDistance: Distance (on axis thickness) from the object to the first surface in mm (float).
        :param primaryWavelength: Primary Wavelength of optical system in mm (float).


        """
        super(OpticalSystem, self).__init__()
        self.surfaces = []
        self.insertSurface(0, Surface( thickness = objectDistance ))  # object
        self.insertSurface(1, Surface())  # image
        # in standard initialization the surface use the BaseAperture which is not limited

        self.primaryWavelength = primaryWavelength


    def appendSurface(self, surface):
        """
        Appends a new surface into the optical system.

        :param position: number of the new surface (int).
           Surface that is currently at this position
           and all following surface indices are incremented.
        """
        self.surfaces.insert(len(self.surfaces)-1, surface)

    def insertSurface(self, position, surface):
        """
        Inserts a new surface into the optical system.

        :param position: number of the new surface (int).
           Surface that is currently at this position
           and all following surface indices are incremented.
        """
        self.surfaces.insert(position, surface)

    def removeSurface(self, position):
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

    def getThickness(self, position):
        """
        Returns the on-axis thickness of a surface.

        :param position: number of the surface (int)
        """
        return self.surfaces[position].getThickness()


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

    def setShape(self, position, shape):
        """
        Sets the shape of a surface.

        :param position: number of the surface (int)
        :param shapeName: name of the Shape child class (str)
        """
        self.surfaces[position].setShape(shape)

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
        #print( abcd )
        return abcd[0, 0] - abcd[0, 1] * abcd[1, 0] / abcd[1, 1]


    def draw2d(self, ax, offset=(0, 0), vertices=100, color="grey"):
        offy = offset[0]
        offz = offset[1]
        for (num, s) in enumerate(self.surfaces):
            s.draw2d(ax, offset=(offy, offz), vertices=vertices, color=color)
            offz += s.getThickness()


    def createOptimizableVariable(self, name, value=0.0, status=False):
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
