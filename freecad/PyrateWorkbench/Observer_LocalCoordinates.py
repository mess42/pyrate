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

import sys
import os
import math

from PySide import QtGui, QtCore
from PySide.QtGui import QInputDialog
from PySide.QtGui import QLineEdit

import FreeCADGui, FreeCAD, Part
from pyrateoptics.raytracer.localcoordinates import LocalCoordinates
from pyrateoptics.core.observers import AbstractObserver

from .Interface_Checks import *

class LC(AbstractObserver):
    def __init__(self, obj, coupling, doc, group):
        if obj == None:
            obj = doc.addObject("Part::FeaturePython", self.returnStructureLabel(coupling.name))

        if group == None:
            group = doc.addObject("App::DocumentObjectGroup", self.returnGroupLabel(coupling.name))

        self.__lc = coupling # link to appropriate data structure
        self.__lc.append_observers([self])
        self.__obj = obj
        self.__group = group
        self.__doc = doc

        group.addObject(obj)
        obj.addProperty("App::PropertyPythonObject", "lcclass", "LC", "class interface").lcclass = coupling
        obj.addProperty("App::PropertyPythonObject", "lcobserver", "LC", "observer interface").lcobserver = self
        obj.addProperty("App::PropertyVector", "globalcoordinates", "LC", "global coords").globalcoordinates = FreeCAD.Base.Vector(tuple(coupling.globalcoordinates))
        obj.addProperty("App::PropertyVector", "decenter", "LC", "decenter").decenter = FreeCAD.Base.Vector((coupling.decx.evaluate(), coupling.decy.evaluate(), coupling.decz.evaluate()))
        obj.addProperty("App::PropertyVector", "tilt", "LC", "tilt in degrees").tilt = FreeCAD.Base.Vector((coupling.tiltx.evaluate()*180.0/math.pi, coupling.tilty.evaluate()*180.0/math.pi, coupling.tiltz.evaluate()*180.0/math.pi))
        obj.addProperty("App::PropertyBool", "order", "LC", "First Tilt then Decenter?").order = bool(coupling.annotations["tiltThenDecenter"])
        obj.addProperty("App::PropertyFloat", "scale", "LC", "Scale factor cross").scale = 1.0

        obj.addProperty("App::PropertyVector", "localbasisX", "LC", "local basis vector x").localbasisX = FreeCAD.Base.Vector(tuple(coupling.localbasis[:,0]))
        obj.addProperty("App::PropertyVector", "localbasisY", "LC", "local basis vector y").localbasisY = FreeCAD.Base.Vector(tuple(coupling.localbasis[:,1]))
        obj.addProperty("App::PropertyVector", "localbasisZ", "LC", "local basis vector z").localbasisZ = FreeCAD.Base.Vector(tuple(coupling.localbasis[:,2]))

        obj.setEditorMode("Placement", 2) # readonly and hide
        obj.setEditorMode("globalcoordinates", 1) # readonly, is determined by tilt, decenter, order
        obj.setEditorMode("localbasisX", 1) # readonly, is determined by tilt, decenter, order
        obj.setEditorMode("localbasisY", 1) # readonly, is determined by tilt, decenter, order
        obj.setEditorMode("localbasisZ", 1) # readonly, is determined by tilt, decenter, order


        obj.Proxy = self

        obj.ViewObject.Proxy=0

        for ch in coupling.getChildren():
            self.createSubgroupForChild(ch)

    def returnGroupLabel(self, s):
        return s + "_LCS"
    def returnStructureLabel(self, s):
        return s

    def createSubgroupForChild(self, lcclasschild):
        subgroup = self.__doc.addObject("App::DocumentObjectGroup", self.returnGroupLabel(lcclasschild.name))
        self.__group.addObject(subgroup)
        chobj = self.__doc.addObject("Part::FeaturePython", self.returnStructureLabel(lcclasschild.name))
        subgroup.addObject(chobj)
        LC(chobj, lcclasschild, self.__doc, subgroup)


    def addChild(self, name=""):
        ch = self.__lc.addChild(LocalCoordinates.p(name))
        self.createSubgroupForChild(ch)


    def getGroup(self):
        return self.__group

    def getLC(self):
        return self.__lc

    group = property(getGroup)

    def inform_about_update(self):
        # override AbstractObserver method
        # let this observer class be informed when update in underlying localcoordinate class takes place
        FreeCAD.Console.PrintMessage("update info from " + self.__lc.name + "\n")
        self.__obj.globalcoordinates = FreeCAD.Base.Vector(tuple(self.__lc.globalcoordinates))
        self.__obj.localbasisX = FreeCAD.Base.Vector(tuple(self.__lc.localbasis[:,0]))
        self.__obj.localbasisY = FreeCAD.Base.Vector(tuple(self.__lc.localbasis[:,1]))
        self.__obj.localbasisZ = FreeCAD.Base.Vector(tuple(self.__lc.localbasis[:,2]))


    def onChanged(self, fp, prop):
        '''Do something when a property has changed'''
        #FreeCAD.Console.PrintMessage("For fp: " + str(fp) + "\n")
        #FreeCAD.Console.PrintMessage("Change property: " + str(prop) + "\n")

        if prop == "Label":
            # rename also the group and the underlying localcoordinate class
            self.__lc.name = self.returnStructureLabel(self.__obj.Label)
            self.__group.Label = self.returnGroupLabel(self.__obj.Label)

        if prop == "order" or prop == "tilt" or prop == "decenter":
            # write back changed properties to underlying localcoordinate class
            # and update tree
            self.__lc.tiltx.setvalue(fp.tilt.x*math.pi/180.0)
            self.__lc.tilty.setvalue(fp.tilt.y*math.pi/180.0)
            self.__lc.tiltz.setvalue(fp.tilt.z*math.pi/180.0)

            self.__lc.decx.setvalue(fp.decenter.x)
            self.__lc.decy.setvalue(fp.decenter.y)
            self.__lc.decz.setvalue(fp.decenter.z)

            self.__lc.tiltThenDecenter = fp.order

            self.__lc.update()

            # perform data structur update of readonly properties after link update
            fp.globalcoordinates = FreeCAD.Base.Vector(tuple(self.__lc.globalcoordinates))
            fp.localbasisX = FreeCAD.Base.Vector(tuple(self.__lc.localbasis[:,0]))
            fp.localbasisY = FreeCAD.Base.Vector(tuple(self.__lc.localbasis[:,1]))
            fp.localbasisZ = FreeCAD.Base.Vector(tuple(self.__lc.localbasis[:,2]))

        if prop != "Shape":
            # prevent recursive call
            self.execute(fp)


    def execute(self, fp):
        '''Do something when doing a recomputation, this method is mandatory'''

        p0 = fp.globalcoordinates
        p1 = p0.add(fp.localbasisX*fp.scale)
        p2 = p0.add(fp.localbasisY*fp.scale)
        p3 = p0.add(fp.localbasisZ*fp.scale)

        l1 = Part.makeLine(p0, p1)
        l2 = Part.makeLine(p0, p2)
        l3 = Part.makeLine(p0, p3)

        c1 = Part.makeCone(0.05*fp.scale, 0.0, 0.1*fp.scale, p1, fp.localbasisX, 360)
        c2 = Part.makeCone(0.05*fp.scale, 0.0, 0.1*fp.scale, p2, fp.localbasisY, 360)
        c3 = Part.makeCone(0.05*fp.scale, 0.0, 0.1*fp.scale, p3, fp.localbasisZ, 360)

        l1.Placement = fp.Placement
        l2.Placement = fp.Placement
        l3.Placement = fp.Placement


        fp.Shape = Part.makeCompound([l1, l2, l3, c1, c2, c3])
        #FreeCAD.Console.PrintMessage("Recompute Python LC feature\n")

    def __getstate__(self):
        '''When saving the document this object gets stored using Python's json module.\
                Since we have some un-serializable parts here -- the Coin stuff -- we must define this method\
                to return a tuple of all serializable objects or None.'''
        return None

    def __setstate__(self,state):
        '''When restoring the serialized object from document we have the chance to set some internals here.\
                Since no data were serialized nothing needs to be done here.'''
        return None



