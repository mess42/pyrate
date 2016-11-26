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



import FreeCAD
import FreeCADGui

from PySide.QtGui import QLineEdit, QInputDialog, QMessageBox


from Observer_OpticalSystem import OpticalSystemObserver 
from Object_MaterialCatalogue import MaterialCatalogueObject
from Object_Material import MaterialObject

from Interface_Identifiers import *
from Interface_Checks import *

from TaskPanel_SurfaceList_Edit import SurfaceListTaskPanelEdit

class CreateSystemTool:
    "Tool for creating optical system"

    def GetResources(self):
        return {"Pixmap"  : ":/icons/pyrate_logo_icon.svg", # resource qrc file needed, and precompile with python-rcc
                "MenuText": "Create optical system ...",
                "Accel": "",
                "ToolTip": "Opens dialog for system creation"
                }

    def IsActive(self):
        if FreeCAD.ActiveDocument == None:
            return False
        else:
            return True

    def Activated(self):

        doc = FreeCAD.ActiveDocument
        (text, ok) = QInputDialog.getText(None, Title_MessageBoxes, "Name for optical system?", QLineEdit.Normal)        
        if text and ok:
            OpticalSystemObserver(doc, text) 

        if existsStandardMaterials(doc):
            result = QMessageBox.question(None, Title_MessageBoxes, 
                                          "No Standard Materials Catalogue defined. Create one?", 
                                          QMessageBox.Yes | QMessageBox.No)
            if result == QMessageBox.Yes:
                stdmatcatalogue = MaterialCatalogueObject(doc, Group_StandardMaterials_Label)
                stdmatcatalogue.addMaterial("Mirror", "Mirror")
                stdmatcatalogue.addMaterial("ConstantIndexGlass", "PMMA", index=1.5)
                stdmatcatalogue.addMaterial("ConstantIndexGlass", "Vacuum")
                stdmatcatalogue.addMaterial("ModelGlass", "mydefaultmodelglass")                

        
        
        # TODO: 1 OSinterface per doc, but several optical systems
        # TODO: call wizard for creation of new system

        # old code

        # doc.OSinterface.dummycreate4() # substitute by system creation dialog
        # dummycreate() -> lens system
        # dummycreate2() -> mirror system
        # dummycreate3() -> lens system with incorrect curvature in surface7
        # dummycreate4() -> GRIN medium
        # doc.OSinterface.createSurfaceViews(doc)
        # doc.OSinterface.showAimFiniteSurfaceStopDialog()
        # doc.OSinterface.showFieldWaveLengthDialog()
        # doc.OSinterface.createRayViews(doc, 50)
        #PyrateInterface.OSinterface.showSpotDiagrams(100)


class ShowSurfaceList:
    def GetResources(self):
        return {"Pixmap"  : ":/icons/pyrate_shape_icon.svg", 
                "MenuText": "Edit Surface List ...",
                "Accel": "",
                "ToolTip": "Manage Surface List"
                }
    
    def IsActive(self):

        if FreeCAD.ActiveDocument == None:
            return False
        else:
            selection = FreeCADGui.Selection.getSelection()
            if len(selection) == 1 and isOpticalSystemObserver(selection[0]): #('wavelengths' in selection[0].PropertiesList):
                # TODO: comparison with CheckObjects function?                
                return True
            else:
                return False

    def Activated(self):
        osselection = FreeCADGui.Selection.getSelection()[0] # only active if len = 1 and obj is appropriate


        # TODO: In menu: Add, Del Surf001, Delf Surf002, ..., Del Surf00N
        # TODO: Enums must have actualized with Labels of defined surfaces
        # TODO: every time something is changed the list in the core opticalsystem has to be actualized
        # TODO: initial Surfaces list from core opticalsystem
        # TODO: initial enumeration list from core opticalsystem
        
        panel = SurfaceListTaskPanelEdit(osselection)
        FreeCADGui.Control.showDialog(panel)



FreeCADGui.addCommand('CreateSystemCommand', CreateSystemTool())
FreeCADGui.addCommand('ShowSurfaceDialogCommand', ShowSurfaceList())

