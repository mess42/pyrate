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
from pyrateoptics.raytracer.surface_shape import (Conic,
                                                  Cylinder,
                                                  Asphere,
                                                  ExplicitShape)
from pyrateoptics.raytracer.aperture import BaseAperture, CircularAperture

import FreeCAD


class SurfaceObject(AbstractObserver):

    def __init__(self, doc, group, name, shapetype, aptype, lclabel,
                 matlabel, **kwargs):
        self.__doc = doc  # e.g. ActiveDocument
        self.__group = group  # surface group
        self.__obj = doc.addObject("Part::FeaturePython", name)
        self.__group.addObject(self.__obj)

        if "surface" in kwargs:
            # initialization from pre defined optical system
            # override shapetype, aptype, lclabel, matlabel could be set later
            surf = kwargs["surface"]
            # self.shape = shape
            # self.material = material
            # self.aperture = aperture
            # self.lc = lc # reference to local coordinate system tree

            shapetype = ""
            if isinstance(surf.shape, Conic):
                shapetype = Shape_Conic
            if isinstance(surf.shape, Cylinder):
                shapetype = Shape_Cylinder
            if isinstance(surf.shape, Asphere):
                shapetype = Shape_Asphere

            aptype = ""
            if isinstance(surf.aperture, BaseAperture):
                aptype = Aperture_Base
            if isinstance(surf.aperture, CircularAperture):
                aptype = Aperture_Circular



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

        self.writebackshapefunc = { # for shape
            Shape_Conic:self.writebackShConic,
            Shape_Cylinder:self.writebackShCylinder,
            Shape_Asphere:self.writebackShAsphere,
            Shape_Explicit:self.writebackShExplicit
        }

        self.readshapefunc = { # for shape
            Shape_Conic:self.readShConic,
            Shape_Cylinder:self.readShCylinder,
            Shape_Asphere:self.readShAsphere,
            Shape_Explicit:self.readShExplicit
        }

        # TODO: set values from initialized matclass coming
        # from a predefined optical system

        self.__obj.addProperty("App::PropertyPythonObject",
                               "shapeclass",
                               "Shape",
                               "surface shape class from pyrateoptics code")

        self.__obj.addProperty("App::PropertyString",
                               "shapetype",
                               "Shape",
                               "specifies type").shapetype = shapetype

        self.__obj.addProperty("App::PropertyPythonObject",
                               "apertureclass",
                               "Aperture",
                               "aperture class from pyrateoptics code")

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

        self.__obj.shapeclass.append_observers([self])
        self.__obj.Proxy = self

    def getObject(self):
        return self.__obj

    # shape initialization

    def initShConic(self, curv=0., cc=0., **kwargs):
        shapeclass = None
        if "surface" in kwargs:
            shapeclass = kwargs["surface"].shape
            curv = shapeclass.curvature.evaluate()
            cc = shapeclass.conic.evaluate()
        else:
            shapeclass = Conic.p(self.__obj.LocalCoordinatesLink.getLC(),
                                 curv=curv, cc=cc)

        self.__obj.addProperty("App::PropertyFloat", "curv", "Shape",
                               "central curvature").curv = curv
        self.__obj.addProperty("App::PropertyFloat", "cc", "Shape",
                               "conic constant").cc = cc
        self.__obj.shapeclass = shapeclass

    def initShCylinder(self, curv=0., cc=0., **kwargs):

        shapeclass = None
        if "surface" in kwargs:
            shapeclass = kwargs["surface"].shape
            curv = shapeclass.curvature.evaluate()
            cc = shapeclass.conic.evaluate()
        else:
            shapeclass = Cylinder.p(curv=curv, cc=cc)

        self.__obj.addProperty("App::PropertyFloat", "curv", "Shape",
                               "central curvature y").curv = curv
        self.__obj.addProperty("App::PropertyFloat", "cc", "Shape",
                               "conic constant y").cc = cc
        self.__obj.shapeclass = shapeclass

    def initShAsphere(self, curv=0., cc=0., asphereparams=[], **kwargs):

        shapeclass = None
        if "surface" in kwargs:
            shapeclass = kwargs["surface"].shape
            curv = shapeclass.curvature.evaluate()
            cc = shapeclass.conic.evaluate()
        else:
            shapeclass = Asphere.p(curv=curv, cc=cc, acoeffs=asphereparams)
        self.__obj.addProperty("App::PropertyFloat", "curv", "Shape",
                               "central curvature").curv = curv
        self.__obj.addProperty("App::PropertyFloat", "cc", "Shape",
                               "conic constant").cc = cc
        self.__obj.addProperty("App::PropertyFloatList", "asphereparams",
                               "Shape", "aspherical corrections").asphereparams = asphereparams
        self.__obj.shapeclass = shapeclass

    def initShExplicit(self, **kwargs):
        # this function offers the flexibility to use function objects
        # with tunable parameters
        # TODO: implement
        self.__obj.shapeclass = Conic.p(curv=0, cc=0)

    # shape readout

    def readShConic(self):
        self.__obj.curv = self.__obj.shapeclass.curvature.evaluate()
        self.__obj.cc = self.__obj.shapeclass.conic.evaluate()

    def readShCylinder(self):
        self.__obj.curv = self.__obj.shapeclass.curvature.evaluate()
        self.__obj.cc = self.__obj.shapeclass.conic.evaluate()

    def readShAsphere(self):
        self.__obj.curv = self.__obj.shapeclass.curvature.evaluate()
        self.__obj.cc = self.__obj.shapeclass.conic.evaluate()
        # TODO: aspheric corrections

    def readShExplicit(self):
        pass

    # shape writeback

    def writebackShConic(self, fp):
        self.__obj.shapeclass.curvature.set_value(fp.curv)
        self.__obj.shapeclass.conic.set_value(fp.cc)

    def writebackShCylinder(self, fp):
        self.__obj.shapeclass.curvature.set_value(fp.curv)
        self.__obj.shapeclass.conic.set_value(fp.cc)

    def writebackShAsphere(self, fp):
        self.__obj.shapeclass.curvature.set_value(fp.curv)
        self.__obj.shapeclass.conic.set_value(fp.cc)

    def writebackShExplicit(self, fp):
        pass

    # aperture initialization

    def initApBase(self, **kwargs):
        self.__obj.apertureclass = BaseAperture.p(self.__obj.LocalCoordinatesLink.Proxy.getLC())

    def initApCircular(self, semidiameter=1.0, tx=0., ty=0., **kwargs):

        apclass = None
        if "surface" in kwargs:
            apclass = kwargs["surface"].aperture
            semidiameter = apclass.semidiameter
            tx = apclass.tx
            ty = apclass.ty
        else:
            apclass = CircularAperture.p(
                self.__obj.LocalCoordinatesLink.Proxy.getLC(),
                semidiameter=semidiameter, tx=tx, ty=ty)

        self.__obj.addProperty("App::PropertyFloat", "semidiameter",
                               "Aperture", "semidiameter").semidiameter = semidiameter
        self.__obj.addProperty("App::PropertyFloat", "tx",
                               "Aperture", "decentration x").tx = tx
        self.__obj.addProperty("App::PropertyFloat", "ty",
                               "Aperture", "decentration y").ty = ty

        self.__obj.apertureclass = apclass

    def initApUserDefined(self, **kwargs):
        self.__obj.apertureclass = BaseAperture.p(
            self.__obj.LocalCoordinatesLink.Proxy.getLC())

    def onChanged(self, fp, prop):
        FreeCAD.Console.PrintMessage("Changed Surface in GUI " + self.__obj.Name + "\n")
        if prop == "Label":
            # what to do if Label is changed?
            pass

        if prop in Surface_GUIChangeableProperties:
            # write back changed properties to underlying material
            self.writebackshapefunc[self.__obj.shapetype](fp)

        # if prop in Aperture_GUIChangeableProperties:
        #     self.writebackaperturefunc[self.__obj.aperturetype](fp)

    def inform_about_update(self):
        # override AbstractObserver method
        FreeCAD.Console.PrintMessage("Changed Surface in Core "
                                     + self.__obj.Name + "\n")

        self.readshapefunc[self.__obj.shapetype]()
