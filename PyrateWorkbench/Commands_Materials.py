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




import FreeCADGui, FreeCAD


from PySide.QtGui import QLineEdit, QInputDialog

from TaskPanel_Materials_Add import MaterialsTaskPanelAdd
from Object_MaterialCatalogue import MaterialCatalogueObject

from Interface_Checks import *
from Interface_Identifiers import *


class CreateMaterialsTool:
    "Tool for creating materials object"

    def GetResources(self):
        return {"Pixmap"  : ":/icons/pyrate_material_icon.svg", # resource qrc file needed, and precompile with python-rcc
                "MenuText": "Create material ...",
                "Accel": "",
                "ToolTip": "Generates material object in document"
                }

    def IsActive(self):
        if FreeCAD.ActiveDocument == None:
            return False
        else:
            return True

    def Activated(self):

        doc = FreeCAD.ActiveDocument


        panel = MaterialsTaskPanelAdd(doc)#, [oso.Label for oso in osobservers])
        FreeCADGui.Control.showDialog(panel)


class CreateMaterialsCatalogueTool:
    "Tool for creating materials object"

    def GetResources(self):
        return {"Pixmap"  : ":/icons/pyrate_material_catalogue_icon.svg", # resource qrc file needed, and precompile with python-rcc
                "MenuText": "Create material catalogue",
                "Accel": "",
                "ToolTip": "Generates material catalogue group in document"
                }

    def IsActive(self):
        if FreeCAD.ActiveDocument == None:
            return False
        else:
            return True

    def Activated(self):

        doc = FreeCAD.ActiveDocument
        
        (text, ok) = QInputDialog.getText(None, Title_MessageBoxes, "Name for material catalogue?", QLineEdit.Normal)        

        if text and ok:
            MaterialCatalogueObject(doc, text)


            
FreeCADGui.addCommand('CreateMaterialsCommand', CreateMaterialsTool())
FreeCADGui.addCommand('CreateMaterialsCatalogueCommand', CreateMaterialsCatalogueTool())

