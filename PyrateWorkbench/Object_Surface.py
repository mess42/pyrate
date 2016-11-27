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

from Interface_Helpers import *
from Interface_Checks import *
from Interface_Identifiers import *

from core.observers import AbstractObserver
from core.surfShape import Conic, Asphere, ExplicitShape
from core.aperture import BaseAperture, CircularAperture

class SurfaceObject(AbstractObserver):


    def __init__(self, doc, group, name, shapetype, aptype, lclabel, matlabel, **kwargs):
        self.__doc = doc # e.g. ActiveDocument
        self.__group = group # surface group
        self.__obj = doc.addObject("App::FeaturePython", name)
        self.__group.addObject(self.__obj)

        self.initshapedict = {
            Shape_Conic:self.initShConic,
            Shape_Cylinder:self.initShCylinder,
            Shape_Asphere:self.initShAsphere,
            Shape_Explicit:self.initShExplicit
            }
            
        self.initaperturedict = {
            Aperture_Base:self.initApBase,
            Aperture_Circular:self.initApCircular,
            Aperture_UserDefined:self.initApUserDefined
        }
#        
#        self.writebackfunc = { # for shape
#            Material_ConstantIndexGlass:self.writebackConstantIndex,
#            Material_ModelGlass:self.writebackModel,
#            Material_GrinMedium:self.writebackGrin
#        }
#        
#        self.readfunc = { # for shape
#            Material_ConstantIndexGlass:self.readConstantIndex,
#            Material_ModelGlass:self.readModel,
#            Material_GrinMedium:self.readGrin
#        }
          
        # TODO: set values from initialized matclass coming from a predefined optical system
        

        self.__obj.addProperty("App::PropertyPythonObject", 
                               "shapeclass", 
                               "Shape", 
                               "surfShape class from core code")
                               
        self.__obj.addProperty("App::PropertyString", 
                               "shapetype", 
                               "Shape", 
                               "specifies type").shapetype = shapetype

        self.__obj.addProperty("App::PropertyPythonObject", 
                               "apertureclass", 
                               "Aperture", 
                               "aperture class from core code")
                               
        self.__obj.addProperty("App::PropertyString", 
                               "aperturetype", 
                               "Aperture", 
                               "specifies type").aperturetype = aptype


        self.__obj.addProperty("App::PropertyString", 
                               "comment", 
                               "Surface", 
                               "comments").comment = ""

        self.__obj.addProperty("App::PropertyLink", 
                               "LocalCoordinatesLink", 
                               "Coordinates", 
                               "local coordinate system").LocalCoordinatesLink = doc.getObjectsByLabel(lclabel)[0]

        self.__obj.addProperty("App::PropertyLink", 
                               "MaterialLink", 
                               "Material", 
                               "material associated with surface").MaterialLink = doc.getObjectsByLabel(matlabel)[0]


        self.initshapedict[shapetype](**kwargs)
        self.initaperturedict[aptype](**kwargs)
        
        self.__obj.shapeclass.appendObservers([self])
        self.__obj.Proxy = self


    # shape initialization    
    
    def initShConic(self, curv=0., cc=0., **kwargs):
        self.__obj.addProperty("App::PropertyFloat", "curv", "Shape", "central curvature").curv = curv
        self.__obj.addProperty("App::PropertyFloat", "cc", "Shape", "conic constant").cc = cc
        self.__obj.shapeclass = Conic(curv=curv, cc=cc)

    def initShCylinder(self, curv=0., cc=0., **kwargs):
        self.__obj.addProperty("App::PropertyFloat", "curv", "Shape", "central curvature y").curv = curv
        self.__obj.addProperty("App::PropertyFloat", "cc", "Shape", "conic constant y").cc = cc
        self.__obj.shapeclass = Cylinder(curv=curv, cc=cc)
    
    def initShAsphere(self, curv=0., cc=0., asphereparams=[], **kwargs):
        self.__obj.addProperty("App::PropertyFloat", "curv", "Shape", "central curvature").curv = curv
        self.__obj.addProperty("App::PropertyFloat", "cc", "Shape", "conic constant").cc = cc
        self.__obj.addProperty("App::PropertyFloatList", "asphereparams", "Shape", "aspherical corrections").asphereparams = asphereparams
        self.__obj.shapeclass = Asphere(curv=curv, cc=cc, acoeffs=asphereparams)

    def initShExplicit(self, **kwargs):
        # this function offers the flexibility to use function objects with tunable parameters
        # TODO: implement
        self.__obj.shapeclass = Conic(curv=0, cc=0)

    # aperture initialization

    def initApBase(self, **kwargs):
        self.__obj.apertureclass = BaseAperture()
        
    def initApCircular(self, semidiameter=1.0, tx=0., ty=0., **kwargs):
        self.__obj.addProperty("App::PropertyFloat", "semidiameter", "Aperture", "semidiameter").semidiameter = semidiameter
        self.__obj.addProperty("App::PropertyFloat", "tx", "Aperture", "decentration x").tx = tx
        self.__obj.addProperty("App::PropertyFloat", "ty", "Aperture", "decentration y").ty = ty
      
        self.__obj.apertureclass = CircularAperture(semidiameter=semidiameter, tx=tx, ty=ty)        
    
    def initApUserDefined(self, **kwargs):
        self.__obj.apertureclass = BaseAperture()