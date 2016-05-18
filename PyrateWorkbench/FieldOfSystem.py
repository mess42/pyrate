#!/usr/bin/python3
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
import Part
import PartGui
import Points


import PyrateInterface
from PySide import QtGui


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

        PyrateInterface.OSinterface.showAimFiniteSurfaceStopDialog()

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
            return True

    def Activated(self):

        PyrateInterface.OSinterface.showFieldWaveLengthDialog()


FreeCADGui.addCommand('ShowAimDialogCommand',ShowAimDialogCommand())
FreeCADGui.addCommand('ShowFieldDialogCommand',ShowFieldDialogCommand())



