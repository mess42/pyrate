#!/usr/bin/env/python
"""
Pyrate - Optical raytracing based on Python

Copyright (C) 2014-2020
               by     Moritz Esslinger moritz.esslinger@web.de
               and    Johannes Hartung j.hartung@gmx.net
               and    Uwe Lippmann  uwe.lippmann@web.de
               and    Thomas Heinze t.heinze@uni-jena.de
               and    others

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

from .Interface_Helpers import *
from .Interface_Checks import *
from .Interface_Identifiers import *

from pyrateoptics.core.observers import AbstractObserver
from pyrateoptics.raytracer.material.material_isotropic import\
    ConstantIndexGlass, ModelGlass
from pyrateoptics.raytracer.material.material_grin import\
    IsotropicGrinMaterial

import numpy as np

import FreeCAD

class MaterialObject(AbstractObserver):


    def __init__(self, doc, group, name, mattype, **kwargs):
        self.__doc = doc # e.g. ActiveDocument
        self.__group = group # material catalogue
        self.__obj = doc.addObject("App::FeaturePython", name)
        self.__group.addObject(self.__obj)

        self.initfundict = {
            Material_ConstantIndexGlass:self.initConstantIndex,
            Material_ModelGlass:self.initModel,
            Material_Mirror:self.initMirror,
            Material_GrinMedium:self.initGrin
            }

        self.writebackfunc = {
            Material_ConstantIndexGlass:self.writebackConstantIndex,
            Material_ModelGlass:self.writebackModel,
            Material_GrinMedium:self.writebackGrin
        }

        self.readfunc = {
            Material_ConstantIndexGlass:self.readConstantIndex,
            Material_ModelGlass:self.readModel,
            Material_GrinMedium:self.readGrin
        }

        # TODO: set values from initialized matclass coming from a predefined optical system

        self.__obj.addProperty("App::PropertyPythonObject", "matclass", "Material", "material class from pyrateoptics code")
        self.__obj.addProperty("App::PropertyString", "comment", "Material", "comments").comment = ""
        self.__obj.addProperty("App::PropertyString", "mattype", "Material", "specifies type").mattype = mattype

        self.initfundict[mattype](**kwargs)
        self.__obj.matclass.append_observers([self])

        self.__obj.Proxy = self
        # TODO: load/save

    # initfunctions

    def initConstantIndex(self, index=1.0):
        self.__obj.addProperty("App::PropertyFloat", "index", "Material", "constant index").index = index
        self.__obj.matclass = ConstantIndexGlass(index)

    def initModel(self, n0=1.49749699179, A=0.0100998734374, B=0.000328623343942):
        self.__obj.addProperty("App::PropertyFloat", "n0", "Material", "Conrady index").n0 = n0
        self.__obj.addProperty("App::PropertyFloat", "A", "Material", "Conrady A").A = A
        self.__obj.addProperty("App::PropertyFloat", "B", "Material", "Conrady B").B = B
        self.__obj.matclass = ModelGlass([n0, A, B])

    def initMirror(self):
        self.__obj.matclass = Mirror()

    def initGrin(
        self,
        fun=lambda x, y, z: np.ones_like(x),
        dfdx=lambda x, y, z: np.zeros_like(x),
        dfdy=lambda x, y, z: np.zeros_like(x),
        dfdz=lambda x, y, z: np.zeros_like(x),
        bndfunction=lambda x, y, z: np.ones_like(x, dtype=bool),
        ds=0.001,
        energyviolation=0.001
        ):
        self.__obj.addProperty("App::PropertyFloat", "ds", "Integration", "Integration step").ds = ds
        self.__obj.addProperty("App::PropertyFloat", "energyviolation", "Integration", "Energy Violation").energyviolation = energyviolation

        self.__obj.matclass = GrinMaterial(fun, dfdx, dfdy, dfdz, ds, energyviolation, bndfunction)


    def writebackConstantIndex(self, fp):
        self.__obj.matclass.n.setvalue(fp.index)

    def writebackModel(self, fp):
        self.__obj.matclass.n0.setvalue(fp.n0)
        self.__obj.matclass.A.setvalue(fp.A)
        self.__obj.matclass.B.setvalue(fp.B)

    def writebackGrin(self, fp):
        pass

    # TODO: GRIN

    def readConstantIndex(self):
        self.__obj.index = self.__obj.matclass.n.evaluate()

    def readModel(self):
        self.__obj.n0 = self.__obj.matclass.n0.evaluate()
        self.__obj.A = self.__obj.matclass.A.evaluate()
        self.__obj.B = self.__obj.matclass.B.evaluate()

    def readGrin(self):
        pass

    def onChanged(self, fp, prop):
        FreeCAD.Console.PrintMessage("Changed Material in GUI " + self.__obj.Name + "\n")
        if prop == "Label":
            # what to do if Label is changed?
            pass

        if prop in Material_GUIChangeableProperties:
            # write back changed properties to underlying material
            self.writebackfunc[self.__obj.mattype](fp)

    def inform_about_update(self):
        # override AbstractObserver method
        FreeCAD.Console.PrintMessage("Changed Material in Core " + self.__obj.Name + "\n")

        self.readfunc[self.__obj.mattype]()

