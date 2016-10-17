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

@author: Johannes Hartung

Implementation of PySide TreeModel/View according to:
http://stackoverflow.com/questions/17278182/qtreeview-with-custom-items

"""
import sys
import os
import math

from PySide import QtGui, QtCore
from PySide.QtGui import QInputDialog
from PySide.QtGui import QLineEdit

import FreeCADGui, FreeCAD, Part
from core.coordinates import LocalCoordinates
from core.observers import AbstractObserver

class LC(AbstractObserver):
    def __init__(self, obj, coupling, doc, group):
        if obj == None:
            obj = doc.addObject("Part::FeaturePython", self.returnStructureLabel(coupling.name))

        if group == None:
            group = doc.addObject("App::DocumentObjectGroup", self.returnGroupLabel(coupling.name))

        self.__lc = coupling # link to appropriate data structure
        self.__lc.observers.append(self)
        self.__obj = obj
        self.__group = group
        self.__doc = doc

        group.addObject(obj)
        obj.addProperty("App::PropertyPythonObject", "lcclass", "LC", "class interface").lcclass = coupling
        obj.addProperty("App::PropertyPythonObject", "lcobserver", "LC", "observer interface").lcobserver = self
        obj.addProperty("App::PropertyVector", "globalcoordinates", "LC", "global coords").globalcoordinates = FreeCAD.Base.Vector(tuple(coupling.globalcoordinates))
        obj.addProperty("App::PropertyVector", "decenter", "LC", "decenter").decenter = FreeCAD.Base.Vector((coupling.decx.evaluate(), coupling.decy.evaluate(), coupling.decz.evaluate()))
        obj.addProperty("App::PropertyVector", "tilt", "LC", "tilt in degrees").tilt = FreeCAD.Base.Vector((coupling.tiltx.evaluate()*180.0/math.pi, coupling.tilty.evaluate()*180.0/math.pi, coupling.tiltz.evaluate()*180.0/math.pi))
        obj.addProperty("App::PropertyBool", "order", "LC", "First Tilt then Decenter?").order = bool(coupling.tiltThenDecenter)
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
        return s + "_group"
    def returnStructureLabel(self, s):
        return s

    def createSubgroupForChild(self, lcclasschild):
        subgroup = self.__doc.addObject("App::DocumentObjectGroup", self.returnGroupLabel(lcclasschild.name))
        self.__group.addObject(subgroup)
        chobj = self.__doc.addObject("Part::FeaturePython", self.returnStructureLabel(lcclasschild.name))
        subgroup.addObject(chobj)                        
        LC(chobj, lcclasschild, self.__doc, subgroup)


    def addChild(self, name=""):
        ch = self.__lc.addChild(LocalCoordinates(name))
        self.createSubgroupForChild(ch)


    def getGroup(self):
        return self.__group

    group = property(getGroup)

    def informUpdate(self):
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
            # TODO: update of all children?
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


class LocalCoordinatesTreeData(LocalCoordinates):
    def __init__(self, name="", **kwargs):
        super(LocalCoordinatesTreeData, self).__init__(name, **kwargs)
        self.icon = None
        self.index = None

    #---------------------------------------------------------------------------
    def position(self):
        position = 0
        if self.parent is not None:
            count = 0
            for child in self.getChildren():
                if child == self:
                    position = count
                    break
                count += 1
        return position

#-------------------------------------------------------------------------------
class LocalCoordinatesTreeModel(QtCore.QAbstractItemModel):

    #---------------------------------------------------------------------------
    def __init__(self, tree):
        super(LocalCoordinatesTreeModel, self).__init__()
        self.__tree = tree
        self.__current = tree

    #---------------------------------------------------------------------------
    def flags(self, index):
        flag = QtCore.Qt.ItemIsEnabled
        if index.isValid():
            flag |= QtCore.Qt.ItemIsSelectable \
                 | QtCore.Qt.ItemIsUserCheckable \
                 | QtCore.Qt.ItemIsEditable \
                 | QtCore.Qt.ItemIsDragEnabled \
                 | QtCore.Qt.ItemIsDropEnabled
        return flag

    #---------------------------------------------------------------------------
    def index(self, row, column, parent=QtCore.QModelIndex()):
        node = QtCore.QModelIndex()
        if parent.isValid():
            nodeS = parent.internalPointer()
            nodeX = nodeS.children[row]
            if column == 0:
                node = self.__createIndex(row, column, nodeX)
            if column == 1:
                # TODO: write down global coordinates
                pass
        else:
            if column == 0:
                node = self.__createIndex(row, column, self.__tree)
            if column == 1:
                pass
        return node

    #---------------------------------------------------------------------------
    def parent(self, index):
        node = QtCore.QModelIndex()
        if index.isValid():
            nodeS = index.internalPointer()
            parent = nodeS.parent
            if parent is not None:
                node = self.__createIndex(parent.position(), 0, parent)
        return node

    #---------------------------------------------------------------------------
    def rowCount(self, index=QtCore.QModelIndex()):
        count = 1
        node = index.internalPointer()
        if node is not None:
            count = len(node.children)
        return count

    #---------------------------------------------------------------------------
    def columnCount(self, index=QtCore.QModelIndex()):
        return 2

    #---------------------------------------------------------------------------
    def data(self, index, role=QtCore.Qt.DisplayRole):
        data = None
        if role == QtCore.Qt.DisplayRole or role == QtCore.Qt.EditRole:
            node = index.internalPointer()
            data = node.name

        if role == QtCore.Qt.ToolTipRole:
            node = index.internalPointer()
            data = "ToolTip " + node.name

        if role == QtCore.Qt.DecorationRole:
            data = QtGui.QIcon("icon.png")
        return data

    #---------------------------------------------------------------------------
    def setData(self, index, value, role=QtCore.Qt.DisplayRole):
        result = True
        if role == QtCore.Qt.EditRole and value != "":
            node = index.internalPointer()
            node.name = value
            result = True
        return result

    #---------------------------------------------------------------------------
    def __createIndex(self, row, column, node):
        if node.index == None:
            index = self.createIndex(row, column, node)
            node.index = index
            icon = QtGui.QIcon("icon.png")
            b = self.setData(index, icon, QtCore.Qt.DecorationRole)
            b = self.setData(index, "ToolTip "+node.name, QtCore.Qt.ToolTipRole)
        return node.index


class LocalCoordinatesTaskPanel:
    def __init__(self):
        # this will create a Qt widget from our ui file
        fn = os.path.join(os.path.dirname(__file__), 'Qt/lcdialog.ui')        
        self.form = FreeCADGui.PySideUic.loadUi(fn)

        data = LocalCoordinatesTreeData(name="test1")
        lc21 = data.addChild(LocalCoordinatesTreeData(name="lc21", decz=40.0))
        lc22 = data.addChild(LocalCoordinatesTreeData(name="lc22", decz=50.0))
        lc23 = data.addChild(LocalCoordinatesTreeData(name="lc23", decz=60.0))
        lc31 = lc22.addChild(LocalCoordinatesTreeData(name="lc31", decz=60.0))
        lc32 = lc22.addChild(LocalCoordinatesTreeData(name="lc32", decz=60.0))
        lc33 = lc23.addChild(LocalCoordinatesTreeData(name="lc33", decz=60.0))
        lc34 = lc23.addChild(LocalCoordinatesTreeData(name="lc34", decz=60.0))
        
        treeModel = LocalCoordinatesTreeModel(data)
        self.form.treeView.setModel(treeModel)#TreeView(treeModel)
        self.form.treeView.setCurrentIndex(treeModel.index(0, 0))

    def accept(self):
        print("pressed ok")
        FreeCADGui.Control.closeDialog()
        
    def reject(self):
        print("pressed cancel")
        FreeCADGui.Control.closeDialog()
        


class CreateLocalCoordinatesTool:
    "Tool for creating local coordinates"

    def GetResources(self):
        return {"Pixmap"  : ":/icons/pyrate_logo_icon.svg", # resource qrc file needed, and precompile with python-rcc
                "MenuText": "Create local coordinates ...",
                "Accel": "",
                "ToolTip": "Opens dialog for local coordinates"
                }

    def IsActive(self):
        if FreeCAD.ActiveDocument == None:
            return False
        else:
            return True

    def Activated(self):

        gad = FreeCAD.ActiveDocument
        if gad == None:
            return

        origin = LocalCoordinates(name="origin")
        lc21 = origin.addChild(LocalCoordinates(name="lc21", decz=40.0))
        lc22 = origin.addChild(LocalCoordinates(name="lc22", tiltx=0.1, decz=50.0))
        lc23 = origin.addChild(LocalCoordinates(name="lc23", decz=60.0))
        lc31 = lc22.addChild(LocalCoordinates(name="lc31", decz=60.0))
        lc32 = lc22.addChild(LocalCoordinates(name="lc32", decz=70.0))
        lc33 = lc23.addChild(LocalCoordinates(decz=60.0))
        lc34 = lc23.addChild(LocalCoordinates(name="lc34", decz=60.0))
    
        llc = LC(None, origin, gad, None)

class ContextAddChildToLocalCoordinatesTool:
    
    "Tool for adding child to local coordinates within context menu"

    def GetResources(self):
        return {"Pixmap"  : ":/icons/pyrate_logo_icon.svg", # resource qrc file needed, and precompile with python-rcc
                "MenuText": "Add child to local coordinates ...",
                "Accel": "",
                "ToolTip": "Add child to local coordinates"
                }


    def IsActive(self):
        if FreeCAD.ActiveDocument == None:
            return False
        else:
            return True

    def Activated(self):
        
        selection = [s  for s in FreeCADGui.Selection.getSelection() if s.Document == FreeCAD.ActiveDocument ]
        (name_of_child, accepted) = QInputDialog.getText(None, "Pyrate", "Name of Child Local Coordinates System", QLineEdit.Normal, "")
        if len(selection) == 1 and accepted:
            obj = selection[0]
            if 'lcclass' in obj.PropertiesList:
                obj.lcobserver.addChild(name = name_of_child)
                
class ContextIncreaseScaleOfAllLocalCoordinatesTool:
    def GetResources(self):
        return {"Pixmap"  : ":/icons/pyrate_logo_icon.svg", # resource qrc file needed, and precompile with python-rcc
                "MenuText": "Increase Scale of local coordinates ...",
                "Accel": "",
                "ToolTip": "increase Scale of local coordinates"
                }

    def IsActive(self):
        if FreeCAD.ActiveDocument == None:
            return False
        else:
            return True

    def Activated(self):

        for o in FreeCAD.ActiveDocument.Objects:
            if "lcclass" in o.PropertiesList:
                o.scale += 1

class ContextDecreaseScaleOfAllLocalCoordinatesTool:
    def GetResources(self):
        return {"Pixmap"  : ":/icons/pyrate_logo_icon.svg", # resource qrc file needed, and precompile with python-rcc
                "MenuText": "Decrease Scale of local coordinates ...",
                "Accel": "",
                "ToolTip": "Decrease Scale of local coordinates"
                }

    def IsActive(self):
        if FreeCAD.ActiveDocument == None:
            return False
        else:
            return True

    def Activated(self):

        for o in FreeCAD.ActiveDocument.Objects:
            if "lcclass" in o.PropertiesList:
                o.scale -= 1
 


FreeCADGui.addCommand('CreateLocalCoordinatesCommand', CreateLocalCoordinatesTool())
FreeCADGui.addCommand('ContextAddChildToLocalCoordinatesCommand', ContextAddChildToLocalCoordinatesTool())
FreeCADGui.addCommand('ContextIncreaseScaleOfAllLocalCoordinatesCommand', ContextIncreaseScaleOfAllLocalCoordinatesTool())
FreeCADGui.addCommand('ContextDecreaseScaleOfAllLocalCoordinatesCommand', ContextDecreaseScaleOfAllLocalCoordinatesTool())