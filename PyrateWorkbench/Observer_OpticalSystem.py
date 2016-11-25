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


# standard include

import time
import math
import sys
import os

import uuid

import numpy as np
import matplotlib.pyplot as plt
from PySide import QtCore, QtGui

from core import material
from core import surfShape
from core import aim
from core import field
from core import pupil
from core import raster
from core import plots
from core import aperture

from core.ray import RayPath
from core.optical_system import OpticalSystem, Surface
from core.observers import AbstractObserver
from core.coordinates import LocalCoordinates
from core.aperture import CircularAperture

# freecad modules

import FreeCAD
import FreeCADGui
import Part
import Points
import Draft


from Observer_LocalCoordinates import LC

from Interface_Identifiers import *
from Interface_Helpers import *
from Interface_Checks import *
    

class OpticalSystemObserver(AbstractObserver):
    def __init__(self, doc, name):
        self.__doc = doc
        obj = doc.addObject("App::FeaturePython", name)
        self.__obj = obj
        obj.Proxy = self

        self.__NameOSGroup = Group_OS_Label + "_" + uuidToName(uuid.uuid4())
        self.__NameSurfaceGroup = Group_Surface_Label + "_" + uuidToName(uuid.uuid4())        
        self.__NameFunctionsGroup = Group_Functions_Label + "_" + uuidToName(uuid.uuid4())
        self.__NameCoordinatesGroup = Group_Coordinates_Label + "_" + uuidToName(uuid.uuid4())
        
        self.__group = doc.addObject("App::DocumentObjectGroup", self.__NameOSGroup)
        self.__group.addObject(obj)
        
        self.__surfacegroup = doc.addObject("App::DocumentObjectGroup", self.__NameSurfaceGroup)
        self.__functionsgroup = doc.addObject("App::DocumentObjectGroup", self.__NameFunctionsGroup)
        self.__coordinatesgroup = doc.addObject("App::DocumentObjectGroup", self.__NameCoordinatesGroup)

        
        self.__group.addObject(self.__surfacegroup)
        self.__group.addObject(self.__functionsgroup)
        self.__group.addObject(self.__coordinatesgroup)

        self.__functionsgroup.Label = Group_Functions_Label + "_" + name
        self.__surfacegroup.Label = Group_Surface_Label + "_" + name
        self.__coordinatesgroup.Label = Group_Coordinates_Label + "_" + name
        self.__group.Label = Group_OS_Label + "_" + name

        
        # TODO: all properties are not really operational

        # group links

        obj.addProperty("App::PropertyString", "NameOSGroup", "Groups", "Name of OS Group").NameOSGroup = self.__NameOSGroup
        obj.addProperty("App::PropertyString", "NameFunctionsGroup", "Groups", "Name of Functions Group").NameFunctionsGroup = self.__NameFunctionsGroup
        obj.addProperty("App::PropertyString", "NameSurfaceGroup", "Groups", "Name of Surface Group").NameSurfaceGroup = self.__NameSurfaceGroup
        obj.addProperty("App::PropertyString", "NameCoordinatesGroup", "Groups", "Name of Coordinates Group").NameCoordinatesGroup = self.__NameCoordinatesGroup
        
        obj.setEditorMode("NameOSGroup", 1) # readonly 
        obj.setEditorMode("NameFunctionsGroup", 1) # readonly
        obj.setEditorMode("NameSurfaceGroup", 1) # readonly
        obj.setEditorMode("NameCoordinatesGroup", 1) # readonly


        # OS Properties
    
        obj.addProperty("App::PropertyPythonObject", 
                        "osclass", 
                        "OS", 
                        "os class interface").osclass = OpticalSystem()

        s = obj.osclass
                        
        lc1 = s.addLocalCoordinateSystem(LocalCoordinates(name="surf1", decz=2.0)) # objectDist
        lc2 = s.addLocalCoordinateSystem(LocalCoordinates(name="surf2", decz=3.0))
        lc3 = s.addLocalCoordinateSystem(LocalCoordinates(name="surf3", decz=5.0, tiltx=0.0*math.pi/180.0))
        lc4 = s.addLocalCoordinateSystem(LocalCoordinates(name="surf4", decz=3.0))
        lc5 = s.addLocalCoordinateSystem(LocalCoordinates(name="surf5", decz=3.0))
        lc6 = s.addLocalCoordinateSystem(LocalCoordinates(name="surf6", decz=2.0))
        lc7 = s.addLocalCoordinateSystem(LocalCoordinates(name="surf7", decz=3.0))
        lc8 = s.addLocalCoordinateSystem(LocalCoordinates(name="image", decz=19.0))
        
        
        s.insertSurface(1, Surface(lc1, surfShape.Conic(curv=1/-5.922), # thickness=3.0,
                                   material=material.ConstantIndexGlass(1.7), 
                                    aperture=CircularAperture(0.55)))
        
        s.insertSurface(2, Surface(lc2, surfShape.Conic(curv=1/-3.160), # thickness=5.0, 
                                   aperture=CircularAperture(1.0)))
        
        s.insertSurface(3, Surface(lc3, surfShape.Conic(curv=1/15.884), #thickness=3.0,
                                   material=material.ConstantIndexGlass(1.7), 
                                    aperture=CircularAperture(1.3)))
        
        s.insertSurface(4, Surface(lc4, surfShape.Conic(curv=1/-12.756), #thickness=3.0,
                                   aperture=CircularAperture(1.3)))
        
        #s.insertSurface(5, Surface(surfShape.Decenter(dx = 0., dy = 1.), material=material.Tilt(angle=20.*np.pi/180.0, axis='X')))
        
        s.insertSurface(5, Surface(lc5, surfShape.Conic(), #thickness=2.0, 
                                   aperture=CircularAperture(1.01))) # Stop Surface
        
        s.insertSurface(6, Surface(lc6, surfShape.Conic(curv=1/3.125), #thickness=3.0,
                                   material=material.ConstantIndexGlass(1.5), 
                                    aperture=CircularAperture(1.0)))
        
        s.insertSurface(7, Surface(lc7, surfShape.Conic(curv=1/1.479), #thickness=19.0,
                                   aperture=CircularAperture(1.0)))
        
        
        s.insertSurface(8, Surface(lc8)) # image
                        
                        
                        
        obj.addProperty("App::PropertyPythonObject", 
                        "coords", 
                        "OS", 
                        "os coords interface").coords = LC(None, obj.osclass.globalcoordinatesystem, doc, self.__coordinatesgroup)
        obj.addProperty("App::PropertyFloatList",
                        "wavelengths",
                        "OS",
                        "wavelengths list").wavelengths = [550.0e-6]

        # Field properties

        obj.addProperty("App::PropertyPythonObject",
                        "fieldpoints",
                        "Field",
                        "Field points").fieldpoints = np.array([[0, 0]])
        obj.addProperty("App::PropertyPythonObject",
                        "fieldpointsbool",
                        "Field",
                        "Field points used?").fieldpointsbool = np.array([True], dtype=bool) 
        obj.addProperty("App::PropertyEnumeration",
                        "fieldtype",
                        "Field",
                        "Type of field?").fieldtype = \
                                ["ObjectHeight",
                                 "ObjectChiefAngle",
                                 "ParaxialImageHeight"]


        # Aiming properties

        obj.addProperty("App::PropertyInteger",
                        "stopposition",
                        "Aiming",
                        "Which surface is stop?").stopposition = 0
        obj.addProperty("App::PropertyEnumeration",
                        "pupiltype",
                        "Aiming",
                        "Type of pupil?").pupiltype = \
                                        ["EntrancePupilDiameter",
                                        "EntrancePupilRadius",
                                        "StopDiameter",
                                        "StopRadius",
                                        "ExitPupilDiameter",
                                        "ExitPupilRadius",
                                        "InfiniteConjugateImageSpaceFNumber",
                                        "InfiniteConjugateObjectSpaceFNumber",
                                        "WorkingImageSpaceFNumber",
                                        "WorkingObjectSpaceFNumber",
                                        "ObjectSpaceNA",
                                        "ImageSpaceNA"]
        obj.addProperty("App::PropertyDistance",
                        "pupilsize",
                        "Aiming",
                        "Pupil size?").pupilsize = 1.0
        obj.addProperty("App::PropertyEnumeration",
                        "rastertype",
                        "Aiming",
                        "Type of pupil rasterization?").rastertype = \
                                                            ["RectGrid",
                                                             "HexGrid",
                                                             "RandomGrid",
                                                             "PoissonDiskSampling",
                                                             "MeridionalFan",
                                                             "SagitalFan",
                                                             "ChiefAndComa",
                                                             "Single"]
                                                             # TODO: -> text file
        obj.addProperty("App::PropertyInteger",
                        "numrays",
                        "Aiming",
                        "How many rays to be drawn?").numrays = 10



    def onChanged(self, fp, prop):
        '''Do something when a property has changed'''
        FreeCAD.Console.PrintMessage("For fp: " + str(fp) + "\n")
        FreeCAD.Console.PrintMessage("Change property: " + str(prop) + "\n")
        
    def informUpdate(self):
        """ can be used if there are any observers coupled to the optical system """
        pass
        
    def __getstate__(self):
        '''When saving the document this object gets stored using Python's json module.\
                Since we have some un-serializable parts here -- the Coin stuff -- we must define this method\
                to return a tuple of all serializable objects or None.'''
        return None
 
    def __setstate__(self,state):
        '''When restoring the serialized object from document we have the chance to set some internals here.\
                Since no data were serialized nothing needs to be done here.'''
        return None

    def getObjectsFromGroupTree(self, grp, boolfun):
        lstboolfun = [o for o in grp.Group if boolfun(o)]
        lstsubgroups = sum([self.getObjectsFromGroupTree(o, boolfun) for o in grp.Group if isGroup(o)], []) # flatten
        return sum([lstboolfun, lstsubgroups], [])
        
        
    
    def returnObjectsFromCoordinatesGroup(self):
        return self.getObjectsFromGroupTree(self.__coordinatesgroup, isLocalCoordinatesObserver)







