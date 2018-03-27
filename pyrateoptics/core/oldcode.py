#!/usr/bin/env/python
"""
Pyrate - Optical raytracing based on Python

Copyright (C) 2017 Moritz Esslinger moritz.esslinger@web.de
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
        return self.material.getABCDMatrix(curvature, nextSurface.getThickness(), nextCurvature, ray)
        # TODO:





class OpticalSystem(ClassWithOptimizableVariables):
    """
    Represents an optical system, consisting of several surfaces and materials inbetween.
    """
    def __init__(self, objectLC = LocalCoordinates(name="object")):
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

    def trace(self, initbundle):
        """
        This function asks the optical system to kindly trace some rays through
        its materials.
        
        :param initbundle (RayBundle object)
        
        :return list of Raybundle objects
        """

        raybundles = [initbundle]

        for (actualSurface, nextSurface) in zip(self.surfaces[:-1], self.surfaces[1:]):
            intersection, t, normal, validIndices = \
                actualSurface.material.propagate(actualSurface, \
                                                nextSurface, \
                                                raybundles[-1])

            raybundles.append(nextSurface.material.refract(actualSurface.material, raybundles[-1], intersection, normal, validIndices))

        return raybundles
        
        
        
        
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
            
            
class ObjectHeight(object):
    def __init__(self):
        pass

    def getChiefSlope(self, opticalSystem, stopPosition, ray, objFieldXY):
        """
        Calculates the chief ray slope from an object field height.

        :param opticalSystem: OpticalSystem object
        :param stopPosition: index of stop surface (int)
        :param ray: raybundle object
        :param objFieldXY: object field height in x and y direction (1d numpy array of 2 floats)

        :return chiefSlopeXY: chief ray slope in x and y direction (1d numpy array of 2 floats)
        """
    
        zen, magen, zex, magex, abcd_obj_stop, abcd_stop_im = opticalSystem.getParaxialPupil(stopPosition, ray)
        chiefSlopeXY = - objFieldXY / zen
        return chiefSlopeXY

    def getObjectHeight(self, opticalSystem, ray, stopPosition, objFieldXY):
        return objFieldXY


class ObjectChiefAngle(ObjectHeight):
    def getChiefSlope(self, opticalSystem, stopPosition, ray, objChiefAngle):
        """
        Calculates the chief ray slope from the object sided chief ray angle.

        :param opticalSystem: OpticalSystem object
        :param ray: raybundle object
        :param objChiefAngle: object sided chief ray angle in degree (1d numpy array of 2 floats)

        :return chiefSlopeXY: chief ray slope in x and y direction (1d numpy array of 2 floats)
        """
    
        return tan(objChiefAngle * pi / 180.0)

    def getObjectHeight(self, opticalSystem, ray, stopPosition, objChiefAngle):
        zen, magen, zex, magex, abcd_obj_stop, abcd_stop_im = opticalSystem.getParaxialPupil(stopPosition, ray)
        objFieldXY = -zen * tan(objChiefAngle * pi / 180.0)
        return objFieldXY


class ParaxialImageHeight(ObjectHeight):
    def getChiefSlope(self, opticalSystem, stopPosition, ray, imFieldXY):
        """
        Calculates the chief ray slope from an image field height assuming no distortion.

        :param opticalSystem: OpticalSystem object
        :param ray: raybundle object
        :param imFieldXY: image field height in x and y direction (1d numpy array of 2 floats)

        :return chiefSlopeXY: chief ray slope in x and y direction (1d numpy array of 2 floats)
        """
    
        zen, magen, zex, magex, abcd_obj_stop, abcd_stop_im = opticalSystem.getParaxialPupil(stopPosition, ray)
        pmag = opticalSystem.getParaxialMagnification(ray)        

        chiefSlopeXY = - imFieldXY / (zen * pmag)
        return chiefSlopeXY

    def getObjectHeight(self, opticalSystem, ray, stopPosition, imFieldXY):
        pmag = opticalSystem.getParaxialMagnification(ray)        
        return imFieldXY / pmag

