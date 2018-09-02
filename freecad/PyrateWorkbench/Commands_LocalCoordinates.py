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

from PySide.QtGui import QInputDialog
from PySide.QtGui import QLineEdit

import FreeCADGui, FreeCAD

from .Interface_Checks import *

from .TaskPanel_LocalCoordinates_Add import LocalCoordinatesTaskPanelAdd


class CreateLocalCoordinatesTool:
    "Tool for creating local coordinates"

    def GetResources(self):
        return {"Pixmap"  : ":/icons/pyrate_coord_icon.svg", # resource qrc file needed, and precompile with python-rcc
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
        if doc == None:
            return

        osobservers = []
        for obj in doc.Objects:
            if isOpticalSystemObserver(obj):
                osobservers.append(obj)
            
            

        lcpanel = LocalCoordinatesTaskPanelAdd(doc, [oso.Label for oso in osobservers])
        FreeCADGui.Control.showDialog(lcpanel)



class ContextAddChildToLocalCoordinatesTool:
    
    "Tool for adding child to local coordinates within context menu"

    def GetResources(self):
        return {"Pixmap"  : ":/icons/pyrate_coord_icon.svg", # resource qrc file needed, and precompile with python-rcc
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
            if isLocalCoordinatesObserver(obj):
                obj.lcobserver.addChild(name = name_of_child)
                
class ContextIncreaseScaleOfAllLocalCoordinatesTool:
    def GetResources(self):
        return {"Pixmap"  : ":/icons/pyrate_coord_icon.svg", # resource qrc file needed, and precompile with python-rcc
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
            if isLocalCoordinatesObserver(o):
                o.scale += 1

class ContextDecreaseScaleOfAllLocalCoordinatesTool:
    def GetResources(self):
        return {"Pixmap"  : ":/icons/pyrate_coord_icon.svg", # resource qrc file needed, and precompile with python-rcc
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
            if isLocalCoordinatesObserver(o):
                o.scale -= 1
 


FreeCADGui.addCommand('CreateLocalCoordinatesCommand', CreateLocalCoordinatesTool())
FreeCADGui.addCommand('ContextAddChildToLocalCoordinatesCommand', ContextAddChildToLocalCoordinatesTool())
FreeCADGui.addCommand('ContextIncreaseScaleOfAllLocalCoordinatesCommand', ContextIncreaseScaleOfAllLocalCoordinatesTool())
FreeCADGui.addCommand('ContextDecreaseScaleOfAllLocalCoordinatesCommand', ContextDecreaseScaleOfAllLocalCoordinatesTool())

