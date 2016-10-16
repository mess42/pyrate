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

import uuid


class Surface(ClassWithOptimizableVariables):
    """
    Represents a surface of an optical system.

    :param shape: Shape of the surface. Calculates the intersection with rays. ( Shape object or child )
    :param material: Material of the volume behind the surface. Calculates the refraction. ( Material object or child )
    :param thickness: distance to next surface on the optical axis
    """
    def __init__(self, lc, shape=surfShape.Conic(), material=ConstantIndexGlass(), aperture=aperture.BaseAperture(), **kwargs):
        super(Surface, self).__init__()

        self.shape = shape
        self.material = material
        self.aperture = aperture
        self.lc = lc # reference to local coordinate system tree

    
    # TODO: these functions will be obsolete, since the thickness parameters is
    # superceded by self.localcoordinates.globalcoordinates and
    # self.localcoordinates.localbasissystem
    def setThickness(self, thickness):
        self.lc.dict_variables["decz"].setvalue(thickness)

    def getThickness(self):
        return self.lc.dict_variables["decz"].evaluate()
        

    def setMaterial(self, material):
        """
        Sets the material object self.mater

        :param material: (object)

        :return self.material: new Material object
        """

        # TODO: conserve most basic material properties

        self.material = material

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
    def __init__(self, objectLC = LocalCoordinates(name="object"), primaryWavelength = 550e-6):
        """
        Creates an optical system object. Initially, it contains 2 plane surfaces (object and image).

        :param objectDistance: Distance (on axis thickness) from the object to the first surface in mm (float).
        :param primaryWavelength: Primary Wavelength of optical system in mm (float).


        """
        super(OpticalSystem, self).__init__()
        
        self.globalcoordinatesystem = LocalCoordinates(name="global")
        self.lcfocus = "global"
        
        self.objectlc = self.addLocalCoordinateSystem(objectLC)
        self.lcfocus = "object"
        

        
        self.surfaces = []
        self.insertSurface(0, Surface(self.objectlc))  # object
        # in standard initialization the surface use the BaseAperture which is not limited

        self.primaryWavelength = primaryWavelength

        #self.observers = {} # observers which will we informed upon change of OS

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
        
    # TODO: removed observer functionality, because observers should be an integrated
    # part of classwithoptimizable variables and maybe at higher levels in the class hierarchy

    #def addObserver(self, name, observer):
    #    self.observers[name] = observer
    
    #def returnObserver(self, name):
    #    return self.observers[name]
        
    #def removeObserver(self, name):
    #    return self.observers.pop(name)
        
    #def informObservers(self):
    #    for o in self.observers:
    #        o.setValues(self.obtainGeometricalSurfaceData())

    #def obtainGeometricalSurfaceData(self):
    #    doublelist = [[s.localcoordinates.thickness.evaluate(), \
    #      s.localcoordinates.decx.evaluate(), \
    #      s.localcoordinates.decy.evaluate(), \
    #      s.localcoordinates.tiltx.evaluate(), \
    #      s.localcoordinates.tilty.evaluate(), \
    #      s.localcoordinates.tiltz.evaluate()] for s in self.surfaces]
    #    return np.array(doublelist)


    def appendSurface(self, surface):
        """
        Appends a new surface into the optical system.

        :param position: number of the new surface (int).
           Surface that is currently at this position
           and all following surface indices are incremented.
        """
        self.surfaces.insert(len(self.surfaces), surface)

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
        # TODO: update of local coordinate references missing
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
        print abcd
        return abcd[0, 0] - abcd[0, 1] * abcd[1, 0] / abcd[1, 1]


    def draw2d(self, ax, vertices=100, color="grey"):
        for (num, s) in enumerate(self.surfaces):
            s.draw2d(ax, vertices=vertices, color=color)
            


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


if __name__ == "__main__":
    os = OpticalSystem()
    lc1 = os.addLocalCoordinateSystem(LocalCoordinates(decz=10.0))
    os.addLocalCoordinateSystem(LocalCoordinates(decz=20.0))
    os.addLocalCoordinateSystem(LocalCoordinates(decz=30.0))
    os.addLocalCoordinateSystem(LocalCoordinates(decz=40.0))
    
    os.addLocalCoordinateSystem(LocalCoordinates(name="COM", decx=10.0, decy=5.0, decz=10.), refname=lc1.name)
    
    print(os.globalcoordinatesystem.pprint())
        
