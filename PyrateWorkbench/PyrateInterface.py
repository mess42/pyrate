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


# freecad modules

import FreeCAD
import FreeCADGui
import Part
import Points
import Draft


from LocalCoordinatesTree import LC


    

class FreeCADOutputStream(object):
    def write(self, txt):
        FreeCAD.Console.PrintMessage(txt)

# TODO: new implementation of OS coupling
class OpticalSystemObserver(AbstractObserver):
    def __init__(self, doc):
        self.__doc = doc
        obj = doc.addObject("App::FeaturePython", "OS")
        self.__obj = obj
        self.__group = doc.addObject("App::DocumentObjectGroup", "OS_group")
        self.__group.addObject(obj)
        self.__surfacegroup = doc.addObject("App::DocumentObjectGroup", "Surfaces_group")
        self.__group.addObject(self.__surfacegroup)

        # TODO: all properties are not really operational

        # OS Properties
    
        obj.addProperty("App::PropertyPythonObject", 
                        "osclass", 
                        "OS", 
                        "os class interface").osclass = OpticalSystem()
        obj.addProperty("App::PropertyPythonObject", 
                        "coords", 
                        "OS", 
                        "os coords interface").coords = LC(None, obj.osclass.globalcoordinatesystem, doc, self.__group)
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








