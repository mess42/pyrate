#!/usr/bin/env/python

import inspect
import shape as surfShape # the name 'shape' already denotes the dimensions of a numpy array
import material

from numpy import *
from ray import RayBundle


class Surface():
    """
    Represents a surface of an optical system.
    
    :param shap: Shape of the surface. Claculates the intersection with rays. ( Shape object or child )
    :param mater: Material of the volume behind the surface. Claculates the refraction. ( Material object or child )
    :param thickness: distance to next surface on the optical axis
    """
    def __init__(self, thickness = 0.0):
        self.shap  = surfShape.Conic()
        self.mater = material.ConstantIndexGlass()
        self.thickness = thickness
        
        self.availableShapeNames, self.availableShapeClasses = self.getAvailableShapes()
        self.availableMaterialTypeNames, self.availableMaterialTypeClasses = self.getAvailableMaterialTypes()


    def getAvailableShapes(self):
        """
        Parses shape.py for class definitions. Returns all but shape.Shape
        
        :return listOfShapeNames: List of strings
        :return listOfShapeClasses: List of Class references
        """
        listOfShapeNames = []
        listOfShapeClasses = []
        for name, cla in inspect.getmembers(surfShape):
            fullname = str(cla).strip()
            if fullname.startswith('<class \'shape.') and fullname != '<class \'shape.Shape\'>':
                listOfShapeNames.append( name )
                listOfShapeClasses.append( cla )
        return listOfShapeNames, listOfShapeClasses

    def getAvailableMaterialTypes(self):
        """
        Parses material.py for class defintions. Returns all but material.Material
        
        :return listOfMaterialTypeNames: List of strings
        :return listOfMaterialTypeClasses: List of Class references
        """
        listOfMaterialTypeNames = []
        listOfMaterialTypeClasses = []
        for name, cla in inspect.getmembers(material):
            fullname = str(cla).strip()
            if fullname.startswith('<class \'material.') and fullname != '<class \'material.Material\'>':
                listOfMaterialTypeNames.append( name )
                listOfMaterialTypeClasses.append( cla )
        return listOfMaterialTypeNames, listOfMaterialTypeClasses


    def setMaterial(self, materialType):
        """
        Sets the material object self.mater
        
        :param materialType: name of the material dispersion formula (str)
        """
        
        if materialType in self.availableMaterialTypeNames:
            i = self.availableMaterialTypeNames.index(materialType)
            self.mater = self.availableMaterialTypeClasses[i]()
        else:
            print 'Warning: material type \'', materialType, '\' not found. setMaterial() aborted.'

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
        """
        if shapeName in self.availableShapeNames:
            i = self.availableShapeNames.index(shapeName)

            # conserve the most basic parameters of the shape
            curv = self.shap.curvature 
            semidiam = self.shap.sdia
         
            self.shap = self.availableShapeClasses[i](curv=curv, semidiam=semidiam)
        else:
            print 'Warning: shape \'', materialType, '\' not found. setShape() aborted.'

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
        return self.mater.getABCDMatrix(curvature, self.thickness, nextCurvature, ray)

class OpticalSystem():
    """
    Represents an optical system, consisting of several surfaces and materials inbetween.
    """
    def __init__(self):
        self.surfaces = []
        self.insertSurface(0) # object
        self.insertSurface(1) # image
        self.stopPosition = None        

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
        self.surfaces[position].thickness = thickness

    def getAvailableShapeNames(self):
        """
        Returns a list of valid Shape child class names (list of str)
        """
        return self.surfaces[0].availableShapeNames

    def getAvailableMaterialTypeNames(self):
        """
        Returns a list of valid Material child class names (list of str)
        """
        return self.surfaces[0].availableMaterialTypeNames

    def setMaterial(self, position, materialType):
        """
        Sets the material of a surface.

        :param position: number of the surface (int)
        :param materialType: name of the Material child class (str)
          See OpticalSystem.getAvailableMaterialTypeNames() for details.
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
          See OpticalSystem.getAvailableShapeNames() for details.
        """
        self.surfaces[position].setShape(shapeName)

    def setStopPosition(self, position):
        """
        Sets one surface as the stop. Don't forget to set its semi-diameter.

        :param position: number of the surface (int) 
        """
        self.stopPosition = position

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

    def getParaxialPupil(self, ray):
        """
        Returns the paraxially calculated pupil positions.

        :param ray: Raybundle object
        :return zen: entrance pupil position from object (float)
        :return magen: entrance pupil magnificaction; entrance pupil diameter per stop diameter (float)
        :return zex: exit pupil position from image (float)
        :return magex: exit pupil magnificaction; exit pupil diameter per stop diameter (float)
        """ 
        abcd = self.getABCDMatrix(ray, 0 , self.stopPosition - 1) # object to stop

        zen  = abcd[0,1] / abcd[0,0] # entrance pupil position from object
        magen = 1.0 / abcd[0,0]     

        abcd = self.getABCDMatrix(ray, self.stopPosition, -1) # stop to image

        zex = - abcd[0,1] / abcd[1,1] # exit pupil position from image
        magex = abcd[0,0] - abcd[0,1] * abcd[1,0] / abcd[1,1]

        return zen, magen, zex, magex

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


    def aimInitialRayBundle(self, pupilType, pupilSize, fieldType, fieldSize, wavelength = 0.55):
        """
        Creates and returns a RayBundle object that aims at the optical system pupil.
        Pupil position is estimated paraxially.
        Aiming into the pupil is non-iterative, which means there is no check 
        whether the ray actually hits the stop at the desired position.
        
        :param pupilType: type of the pupil definition (str)
        :param pupilSize: pupil or stop size. Definition of size and unit depends on pupilType, see below. (float)
        :param fieldType: type of the field definition (str)
        :param fieldSize: Field position or angle. Definition and unit depends on fieldType, see below. (float)
        :param wavelength: wavelength in um (float)
        :return raybundle: RayBundle object
       
          pupilType                    definition of pupilSize
          "stop radius"                stop radius in lens units
          "stop diameter"              stop diameter in lens units
          "entrance pupil radius"      pupil radius in lens units
          "entrance pupil diameter"    pupil diameter in lens units
          "exit pupil radius"          pupil radius in lens units
                                         uses paraxial estimate for ratio of entrance and exit pupil size
          "exit pupil diameter"        pupil diameter in lens units
                                         uses paraxial estimate for ratio of entrance and exit pupil size
          "object space f#"            aperture number, unitless
                                         infinite conjugate object sided aperture number
                                         (ratio of paraxial exit pupil diameter and paraxial focal length)
          "image space f#"             aperture number, unitless
                                         infinite conjugate image sided aperture number
                                         (ratio of entrance pupil diameter and paraxial focal length)
          "object space working f#"    aperture number, unitless
                                         finite conjugate object sided aperture number
          "image space working f#"     aperture number, unitless
                                         finite conjugate image sided aperture number
          "object space na"            numerical aperture, unitless
          "image space na"             numerical aperture, unitless
                                         uses paraxial magnification to convert image to object sided NA
    
          fieldType                            definition of fieldSize
          "object height"                      object height in lens units
                                                 Distance from object to first lens must be finite.
          "image height"                       paraxial image height in lens units
                                                 System must not be image sided afocal.
          "infinite conjugate image height"    paraxial image height for incident collimated beam in lens units
                                                 System must not be image sided afocal.
                                                 Thickness of object may be finite.
          "chief ray angle"                    finite conjugate chief ray angle from object to entrance pupil in deg
                                                 System must not be object sided telecentric.
          "infinite conjugate chief ray angle" chief ray angle of collimated incident beam in deg
                                                 Thickness of object may be finite.

        TO DO: complete implementation of pupil definitions
        TO DO: At the moment, this function fails to produce correct values for immersion
        """
        pupilType = pupilType.upper().strip()

        raybundle = RayBundle(zeros((3,3)), zeros((3,3)), wavelength) # dummy ray
        f = self.getEffectiveFocalLength(raybundle)
        zen, magen, zex, magex = self.getParaxialPupil(raybundle)
        pmag = self.getParaxialMagnification(raybundle)        

        # definition of entrance pupil radius epr
        if pupilType == "STOP RADIUS":
            epr = magen * pupilSize
        if pupilType == "STOP RADIUS":
            epr = magen * pupilSize
        if pupilType == "ENTRANCE PUPIL RADIUS":
            epr = pupilSize
        if pupilType == "ENTRANCE PUPIL DIAMETER":
            epr = 0.5 * pupilSize
        if pupilType == "EXIT PUPIL RADIUS":
            epr =  magen / magex * pupilSize
        if pupilType == "EXIT PUPIL DIAMETER":
            epr =  0.5 * magen / magex * pupilSize
        if pupilType == "IMAGE SPACE F#":
            epr = 0.5 * f / pupilSize
        if pupilType == "OBJECT SPACE F#":
            expr = 0.5 * f / pupilSize
            epr =  magen / magex * expr
        if pupilType == "IMAGE SPACE WORKING F#":
            raise NotImplementedError()
        if pupilType == "OBJECT SPACE WORKING F#":
            raise NotImplementedError()
        if pupilType == "OBJECT SPACE NA":
            epr = zen * tan( arcsin( pupilSize ) )
        if pupilType == "IMAGE SPACE NA":
            object_space_na = size * pmag
            epr = zen * tan( arcsin( object_space_na ) )
            # this is only correct for small NA s
            raise NotImplementedError()
        
        return raybundle


    def draw2d(self, ax, offset = [0,0], vertices=100, color="grey"):
        N = self.getNumberOfSurfaces()
        offy = offset[0]
        offz = offset[1]
        for i in arange(N):
            self.surfaces[i].draw2d(ax, offset = [offy, offz])
            offz += self.surfaces[i].thickness
 

