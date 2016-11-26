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
import FreeCADGui

from Interface_Helpers import *
from Interface_Checks import *
from Interface_Identifiers import *

from core.material import ConstantIndexGlass, ModelGlass, Mirror

class MaterialObject:


    def __init__(self, doc, group, name, mattype, **kwargs):
        self.__doc = doc # e.g. ActiveDocument
        self.__group = group # functions group
        self.__obj = doc.addObject("App::FeaturePython", name)
        self.__group.addObject(self.__obj)

        self.initfundict = {
            Material_ConstantIndexGlass:self.initConstantIndex,
            Material_ModelGlass:self.initModel,
            Material_Mirror:self.initMirror,
            }
          
        # TODO: which functions should be called if onChanged event occurs (to write changed values back)
        # TODO: set values from initialized matclass coming from a predefined optical system
        
        #self.changefundict = {
        #}

        self.__obj.addProperty("App::PropertyPythonObject", "matclass", "Material", "material class from core code")
        self.__obj.addProperty("App::PropertyString", "comment", "Material", "comments").comment = ""

        self.initfundict[mattype](**kwargs)

        self.__obj.Proxy = self
        # TODO: load/save

    # initfunctions

    def initConstantIndex(self, index=1.0):
        self.__obj.addProperty("App::PropertyFloat", "index", "Material", "constant index").index = index
        self.__obj.comment = Material_ConstantIndexGlass
        self.__obj.matclass = ConstantIndexGlass(index)

    def initModel(self, n0=1.49749699179, a=0.0100998734374, b=0.000328623343942):
        self.__obj.addProperty("App::PropertyFloat", "n0", "Material", "constant index").n0 = n0
        self.__obj.addProperty("App::PropertyFloat", "a", "Material", "constant index").a = a
        self.__obj.addProperty("App::PropertyFloat", "b", "Material", "constant index").b = b        
        self.__obj.comment = Material_ModelGlass        
        self.__obj.matclass = ModelGlass([n0, a, b])
        
    def initMirror(self):
        self.__obj.comment = Material_Mirror        
        self.__obj.matclass = Mirror()
        
    # TODO: GRIN    
