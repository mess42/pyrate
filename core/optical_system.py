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

import numpy as np

from material import ConstantIndexGlass
from localcoordinates import LocalCoordinates
from localcoordinatestreebase import LocalCoordinatesTreeBase
from ray import RayPath

class OpticalSystem(LocalCoordinatesTreeBase):
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

        super(OpticalSystem, self).__init__(self.rootcoordinatesystem, label = name)
        
        if matbackground is None:
            matbackground = ConstantIndexGlass(self.rootcoordinatesystem, 1.0)

        self.material_background = matbackground # Background material        
        self.elements = {}

    def seqtrace(self, initialbundle, elementsequence): # [("elem1", [1, 3, 4]), ("elem2", [1,4,4]), ("elem1", [4, 3, 1])]
        rpath = RayPath(initialbundle)
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

   
    os = OpticalSystem()
    
    lc1 = os.addLocalCoordinateSystem(LocalCoordinates(decz=10.0), refname=os.rootcoordinatesystem.name)
    lc2 = os.addLocalCoordinateSystem(LocalCoordinates(decz=20.0), refname=lc1.name)
    lc3 = os.addLocalCoordinateSystem(LocalCoordinates(decz=30.0), refname=lc2.name)
    lc4 = os.addLocalCoordinateSystem(LocalCoordinates(decz=40.0), refname=lc3.name)
    
    lc5 = os.addLocalCoordinateSystem(LocalCoordinates(name="COM", decx=10.0, decy=5.0, decz=10.), refname=lc1.name)
    
    print(os.rootcoordinatesystem.pprint())



       
