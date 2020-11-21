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

import FreeCAD
import FreeCADGui

from PySide import QtGui

from .Interface_Identifiers import Title_MessageBoxes
from .TaskPanel_Functions_Edit import FunctionsTaskPanelEdit
from .Object_NotSerializable import NotSerializable

from pyrateoptics.core.functionobject import FunctionObject


class FunctionsView(NotSerializable):

    def __init__(self, vobj):
        self.vobj = vobj
        self.vobj.Proxy = self

        self.obj = vobj.Object
        self.path = ""

    def doubleClicked(self, obj):
        panel = FunctionsTaskPanelEdit(self.obj)
        FreeCADGui.Control.showDialog(panel)

    def setupContextMenu(self, obj, menu):
        actload = menu.addAction("Load source code")
        actsave = menu.addAction("Save source code")
        acttrust = menu.addAction("Trust")
        actuntrust = menu.addAction("Untrust")

        actload.triggered.connect(self.loadFile)
        actsave.triggered.connect(self.saveFile)
        acttrust.triggered.connect(self.trust)
        actuntrust.triggered.connect(self.untrust)

    def trust(self):
        self.obj.trusted = True
        #self.obj.Proxy.trust()

    def untrust(self):
        self.obj.trusted = False
        #self.obj.Proxy.untrust()

    def loadFile(self):
        (FileName, Result) = QtGui.QFileDialog.getOpenFileName(
            None, Title_MessageBoxes, "",
            "Source files (*.py *.FCMacro);;All files (*.*)")

        if Result:
            # sets automatically source_checked flag to false
            self.obj.Proxy.CoreFunctionObject.load(FileName)
            result = QtGui.QMessageBox.question(
                None, Title_MessageBoxes,
                "Did you verify the source code?",
                QtGui.QMessageBox.Yes | QtGui.QMessageBox.No)
            if result == QtGui.QMessageBox.Yes:
                self.trust()
            else:
                self.untrust()

    def saveFile(self):
        QtGui.QMessageBox.information(
            None, Title_MessageBoxes, self.obj.Proxy.CoreFunctionObject.source)
        (FileName, Result) = QtGui.QFileDialog.getSaveFileName(
            None, Title_MessageBoxes, self.path,
            "Source files (*.py *.FCMacro);;All files (*.*)")
        if Result:
            self.obj.Proxy.CoreFunctionObject.save(FileName)

# TODO: adding objects to document and document_group_path
# document_group_path: "/group1/subgroup_a/object_label",
# "/" representing the root in the tree
# "group1" representing a group as well as "subgroup_a"
# "object_label" is the label of the object
# FunctionObjects can only be added to FunctionObjectPools
# In FreeCAD these should be subgroups in the document with certain
# properties


class FunctionsObject(NotSerializable):

    def __init__(self, name, initialsrc, doc, group):
        self.Document = doc  # e.g. ActiveDocument
        self.Group = group  # functions group
        mylabel = self.returnStructureLabel(name)
        self.CoreFunctionObject = FunctionObject(initialsrc, name=mylabel)

        self.Object = doc.addObject("App::FeaturePython", mylabel)
        self.Group.addObject(self.Object)
        self.Object.addProperty("App::PropertyStringList", "functions",
                                "FunctionObject",
                                "functions in object").functions = []
        self.Object.addProperty("App::PropertyBool", "trusted",
                                "FunctionObject",
                                "trusted?").trusted = True

        self.Object.Proxy = self

    def returnFunctionObjects(self):
        self.CoreFunctionObject.generate_functions_from_source(
            self.Object.functions)
        return self.CoreFunctionObject.functions

    def returnSingleFunctionObject(self, name):
        self.CoreFunctionObject.generate_functions_from_source([name])
        result = self.CoreFunctionObject.functions.get(name, None)
        # if dictionary is not filled (due to security issues), return None

        return result

    def trust(self):
        self.CoreFunctionObject.sourcecode_security_checked = True
        self.CoreFunctionObject.globals_security_checked = True

    def untrust(self):
        self.CoreFunctionObject.sourcecode_security_checked = False
        self.CoreFunctionObject.globals_security_checked = False

    def returnStructureLabel(self, name):
        return "function_" + name

    def onChanged(self, fp, prop):
        if prop == "trusted":
            if self.Object.trusted:
                self.trust()
            else:
                self.untrust()
