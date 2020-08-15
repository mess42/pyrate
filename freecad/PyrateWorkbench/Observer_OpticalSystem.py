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

import uuid
import math
from ast import literal_eval

import numpy as np

from pyrateoptics.core.observers import AbstractObserver

from pyrateoptics.raytracer.material import material_isotropic
from pyrateoptics.raytracer.surface_shape import Conic

from pyrateoptics.raytracer.optical_system import OpticalSystem
from pyrateoptics.raytracer.optical_element import OpticalElement
from pyrateoptics.raytracer.surface import Surface
from pyrateoptics.raytracer.localcoordinates import LocalCoordinates

from pyrateoptics.raytracer import ray

# freecad modules

import FreeCAD, Part

# TODO: rename Observer to object (developer sees if it is derived from Observer)
from .Observer_LocalCoordinates import LC

from .Object_Surface import SurfaceObject
from .View_Surface import SurfaceView

from .Interface_Identifiers import *
from .Interface_Helpers import *
from .Interface_Checks import *


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
                        "os class interface").osclass = OpticalSystem.p()






        obj.addProperty("App::PropertyPythonObject",
                        "coords",
                        "OS",
                        "os coords interface").coords = LC(None, obj.osclass.rootcoordinatesystem, doc, self.__coordinatesgroup)
        obj.addProperty("App::PropertyFloatList",
                        "wavelengths",
                        "OS",
                        "wavelengths list").wavelengths = [550.0e-6]
        obj.addProperty("App::PropertyLinkList", "surfaces", "OS", "surface list").surfaces = []


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


    def initDemoSystem(self):
        s = OpticalSystem.p()

        lc0 = s.addLocalCoordinateSystem(LocalCoordinates.p(name="object", decz=0.0), refname=s.rootcoordinatesystem.name)
        lc1 = s.addLocalCoordinateSystem(LocalCoordinates.p(name="surf1", decz=2.0), refname=lc0.name) # objectDist
        lc2 = s.addLocalCoordinateSystem(LocalCoordinates.p(name="surf2", decz=3.0), refname=lc1.name)
        lc3 = s.addLocalCoordinateSystem(LocalCoordinates.p(name="surf3", decz=5.0, tiltx=0.0*math.pi/180.0), refname=lc2.name)
        lc4 = s.addLocalCoordinateSystem(LocalCoordinates.p(name="surf4", decz=3.0), refname=lc3.name)
        lc5 = s.addLocalCoordinateSystem(LocalCoordinates.p(name="surf5", decz=3.0), refname=lc4.name)
        lc6 = s.addLocalCoordinateSystem(LocalCoordinates.p(name="surf6", decz=2.0), refname=lc5.name)
        lc7 = s.addLocalCoordinateSystem(LocalCoordinates.p(name="surf7", decz=3.0), refname=lc6.name)
        lc8 = s.addLocalCoordinateSystem(LocalCoordinates.p(name="image", decz=19.0), refname=lc7.name)


        objectsurf = Surface.p(lc0)
        surf1 = Surface.p(lc1, shape=Conic.p(lc1, curv=1/-5.922))
        surf2 = Surface.p(lc2, shape=Conic.p(lc2, curv=1/-3.160))
        surf3 = Surface.p(lc3, shape=Conic.p(lc3, curv=1/15.884))
        surf4 = Surface.p(lc4, shape=Conic.p(lc4, curv=1/-12.756))
        stopsurf = Surface.p(lc5)
        surf6 = Surface.p(lc6, shape=Conic.p(lc6, curv=1/3.125))
        surf7 = Surface.p(lc7, shape=Conic.p(lc7, curv=1/1.479))
        image = Surface.p(lc8)


        elem = OpticalElement.p(lc0, name="lenssystem")

        glass = material_isotropic.ConstantIndexGlass.p(lc0, n=1.7)
        glass2 = material_isotropic.ConstantIndexGlass.p(lc0, n=1.5)

        elem.addMaterial("glass", glass)
        elem.addMaterial("glass2", glass2)

        elem.addSurface("object", objectsurf, (None, None))
        elem.addSurface("surf1", surf1, (None, "glass"))
        elem.addSurface("surf2", surf2, ("glass", None))
        elem.addSurface("surf3", surf3, (None, "glass"))
        elem.addSurface("surf4", surf4, ("glass", None))
        elem.addSurface("stop", stopsurf, (None, None))
        elem.addSurface("surf6", surf6, (None, "glass2"))
        elem.addSurface("surf7", surf7, ("glass2", None))
        elem.addSurface("image", image, (None, None))

        for mysurf in elem.surfaces.values():
            print(mysurf.aperture.annotations)

        s.addElement("lenssys", elem)

        return s


    def initFromGivenOpticalSystem(self, s):

        # delete surfaces and coordinate systems before fill them up from s
        # do not remove functions objects

        self.__surfacegroup.removeObjectsFromDocument()
        self.__coordinatesgroup.removeObjectsFromDocument()

        # TODO: reference error induced because reference to variables in objects vanishes
        # due to reference counting

        self.__obj.osclass = s
        self.__obj.coords = LC(None, s.rootcoordinatesystem, self.__doc, self.__coordinatesgroup)

        # first init coordinate systems then surfaces

        for (key_elem, elem) in list(s.elements.items()):
            for (key_surf, surf) in list(elem.surfaces.items()):
                so = SurfaceObject(self.__doc, self.__surfacegroup, key_surf, shapetype="", aptype="", lclabel="global", matlabel="Vacuum", surface=surf)
                SurfaceView(so.getObject().ViewObject)
                so.getObject().LocalCoordinatesLink = self.__doc.getObject(surf.rootcoordinatesystem.name) # update local coordinates links


    def makeRayBundle(self, raybundle, offset):
        raysorigin = raybundle.o
        nrays = np.shape(raysorigin)[1]

        pp = Points.Points()
        sectionpoints = []

        res = []

        for i in range(nrays):
            if abs(raybundle.t[i]) > 1e-6:
                x1 = raysorigin[0, i] + offset[0]
                y1 = raysorigin[1, i] + offset[1]
                z1 = raysorigin[2, i] + offset[2]

                x2 = x1 + raybundle.t[i] * raybundle.rayDir[0, i]
                y2 = y1 + raybundle.t[i] * raybundle.rayDir[1, i]
                z2 = z1 + raybundle.t[i] * raybundle.rayDir[2, i]

                res.append(Part.makeLine((x1,y1,z1),(x2,y2,z2))) # draw ray
                sectionpoints.append((x2,y2,z2))
        pp.addPoints(sectionpoints)
        #Points.show(pp) # draw intersection points per raybundle per field point

        return (pp, res)


    def makeRaysFromRayPath(self, raypath, offset, color = (0.5, 0.5, 0.5)):
        def shift(l, n):
            return l[n:] + l[:n]

        #grincolor = (np.random.random(), np.random.random(), np.random.random())

        doc = FreeCAD.ActiveDocument # in initialisierung auslagern
        Nraybundles = len(raypath.raybundles)
        offx = offset[0]
        offy = offset[1]
        offz = offset[2]

        for i in np.arange(Nraybundles):
            offz += self.os.surfaces[i].getThickness()

            (intersectionpts, rays) = self.makeRayBundle(raypath.raybundles[i], offset=(offx, offy, offz))

            #FreeCAD.Console.PrintMessage(str(intersectionpts)+'\n')

            FCptsobj = doc.addObject("Points::Feature", "Surf_"+str(i)+"_Intersectionpoints")
            FCptsobj.Points = intersectionpts
            FCptsview = FCptsobj.ViewObject
            FCptsview.PointSize = 5.0
            FCptsview.ShapeColor = (1.0, 1.0, 0.0)

            self.intersectptsobs.append(FCptsobj)

            for (n, ray) in enumerate(rays):
                FCrayobj = doc.addObject("Part::Feature", "Surf_"+str(i)+"_Ray_"+str(n))
                FCrayobj.Shape = ray
                FCrayview = FCrayobj.ViewObject

                FCrayview.LineColor = color
                FCrayview.PointColor = (1.0, 1.0, 0.0)

                self.rayobs.append(FCrayobj)



    def calculateAimys(self):


        pupiltype = literal_eval("pupil." + self.__obj.pupiltype) # eval is evil but who cares :p
        pupilsize = self.__obj.pupilsize.Value
        fieldType = literal_eval("field."+ self.__obj.fieldtype)
        rasterType = literal_eval("raster."+ self.__obj.rastertype)
        stopPosition = self.__obj.stopposition
        numrays = self.__obj.numrays

        aimys = [aim.aimFiniteByMakingASurfaceTheStop(self.__obj.osclass, pupilType=pupiltype, \
                                                    pupilSizeParameter=pupilsize, \
                                                    fieldType=fieldType, \
                                                    rasterType=rasterType, \
                                                    nray=numrays, wavelength=w, \
                                                    stopPosition=stopPosition) for w in self.__obj.wavelengths]

        return aimys

    def calculateRaypaths(self, aimys):
        raypaths = []
        fieldpoints = self.__obj.fieldpoints[self.__obj.fieldpointsbool]
        for (w, aimy) in zip(self.__obj.wavelengths, aimys):
            for fXY in fieldpoints:
                initialBundle = aimy.getInitialRayBundle(self.__obj.osclass, fieldXY=fXY, wavelength=w) # why we need the wavelength another time?
                raypath = ray.RayPath(initialBundle, self.__obj.osclass)
                raypaths.append(raypath)
        return raypaths

    def drawRaypaths(self, raypaths):
        for rp in raypaths: # per field point
            compoundlist = []
            for (r1, r2) in zip(rp.raybundles[:-1], rp.raybundles[1:]): # per surface
                for (p1, p2) in zip(r1.o.T.tolist(), r2.o.T.tolist()):
                    compoundlist.append(Part.makeLine(
                        FreeCAD.Base.Vector(*p1), FreeCAD.Base.Vector(*p2)
                    ))

            cp = Part.makeCompound(compoundlist)
            FCrayobj = self.__doc.addObject("Part::Feature", "Ray")
            FCrayobj.Shape = cp
            FCrayview = FCrayobj.ViewObject

            FCrayview.LineColor = tuple(np.random.random(3).tolist())
            FCrayview.PointColor = (1.0, 1.0, 0.0)


            #Part.show(cp)

