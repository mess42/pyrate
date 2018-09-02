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

import FreeCAD, FreeCADGui

from .TaskPanel_FieldPoints import FieldPointsTaskPanel
from .Interface_Checks import *


class ShowAimDialogCommand:
    "Show aim dialog"

    def GetResources(self):
        return {"MenuText": "Aim Dialog",
                "Accel": "",
                "ToolTip": "You may enter data for the aiming",
                "Pixmap": ":/icons/pyrate_del_sys_icon.svg"
                }

    def IsActive(self):
        if FreeCAD.ActiveDocument == None:
            return False
        else:
            return True

    def Activated(self):
        pass
        #PyrateInterface.OSinterface.showAimFiniteSurfaceStopDialog()

class ShowFieldDialogCommand:
    "Show field dialog"

    def GetResources(self):
        return {"MenuText": "Field Dialog",
                "Accel": "",
                "ToolTip": "You may enter field point data",
                "Pixmap": ":/icons/pyrate_del_sys_icon.svg"
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
        
        fppanel = FieldPointsTaskPanel(osselection)
        FreeCADGui.Control.showDialog(fppanel)


FreeCADGui.addCommand('ShowAimDialogCommand',ShowAimDialogCommand())
FreeCADGui.addCommand('ShowFieldDialogCommand',ShowFieldDialogCommand())
