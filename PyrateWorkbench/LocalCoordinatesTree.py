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
from PySide import QtGui, QtCore
import FreeCADGui, FreeCAD
from core.coordinates import LocalCoordinates




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

        doc = FreeCAD.ActiveDocument
        panel = LocalCoordinatesTaskPanel()
        FreeCADGui.Control.showDialog(panel)



FreeCADGui.addCommand('CreateLocalCoordinatesCommand', CreateLocalCoordinatesTool())
