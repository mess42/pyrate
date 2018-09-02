#!/usr/bin/env/python
"""
Pyrate - Optical raytracing based on Python

Copyright (C) 2014-2018
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

import FreeCADGui

from PySide import QtGui

from .Interface_Identifiers import *
from .TaskPanel_Functions_Edit import FunctionsTaskPanelEdit


class FunctionsView:
    
    def __init__(self, vobj):
        self.vobj = vobj
        self.vobj.Proxy = self
        
        self.obj = vobj.Object
        self.path = ""

    
    def doubleClicked(self, obj):
        # TODO: open text editor window (dlg_functionobjects_edit.ui)
        # TODO: implement widget class for text editor with little syntax coloring
        # TODO: implement widget class with line numbering
        panel = FunctionsTaskPanelEdit(self.obj)
        FreeCADGui.Control.showDialog(panel)
        
        
        
    def setupContextMenu(self, obj, menu):
        actload = menu.addAction("Load function object")
        actsave = menu.addAction("Save function object")
        
        actload.triggered.connect(self.loadFile)
        actsave.triggered.connect(self.saveFile)        

    def loadFile(self):
        (FileName, Result) = QtGui.QFileDialog.getOpenFileName(None, Title_MessageBoxes, "", "Source files (*.py *.FCMacro);;All files (*.*)")
        
        if Result:
            self.path = FileName
            fp = open(FileName)
            self.obj.Proxy.Source = "".join([line for line in fp])
            fp.close()
    
    def saveFile(self):
        QtGui.QMessageBox.information(None, Title_MessageBoxes,self.obj.Proxy.Source)
        (FileName, Result) = QtGui.QFileDialog.getSaveFileName(None, Title_MessageBoxes, self.path, "Source files (*.py *.FCMacro);;All files (*.*)")        
        # TODO: implement file save procedure, but let messagebox show the first few lines of the code
        if Result:
            fp = open(FileName, "w")
            fp.write(self.obj.Proxy.Source)
            fp.close()

    def __setstate__(self, state):
        return None
        
    def __getstate__(self):
        return None



class FunctionsObject:
    
    
    def __init__(self, name, initialsrc, doc, group):
        self.Document = doc # e.g. ActiveDocument
        self.Group = group # functions group
        self.Object = doc.addObject("App::FeaturePython", self.returnStructureLabel(name))
        self.Group.addObject(self.Object)
        self.Object.addProperty("App::PropertyStringList", "functions", "FunctionObject", "functions in object").functions = []

        self.Object.Proxy = self
        self.Source = initialsrc
        
        # TODO: load/save
        
    def getFunctionsFromSource(self, sourcecodestring, funcnamelist):
        localsdict = {}
        functionsobjects = []
        try:
            exec(sourcecodestring, localsdict) 
            # exec is a security risk, but the code to be executed is loaded from a file at most
            # which has to be inspected by the user; there is no automatic code execution
            
            # TODO: substitute this code by execfile interface
            
        except:
            QtGui.QMessageBox.information(None, Title_MessageBoxes,"Exception caught. Problem in " + self.Object.Label)
            return functionsobjects
            
        for fn in funcnamelist:
            try:
                functionsobjects.append(localsdict[fn])
            except:
                pass
        return functionsobjects
        
    def createSourceCode(self):
        return str(self.Source) 
        # removed reference to stringlist property source, due to security risk
        # of loading and saving code from/to FreeCAD documents
    
    
    def returnFunctionObjects(self):
        return self.getFunctionsFromSource(self.createSourceCode(), self.Object.functions)
        
    def returnSingleFunctionObject(self, name):
        result = None
        try:
            result = self.getFunctionsFromSource(self.createSourceCode(), [name])[0]
        except:
            pass
        
        return result
        
        
    def returnStructureLabel(self, name):
        return "function_" + name

    
    def __setstate__(self, state):
        return None
        
    def __getstate__(self):
        return None
        
