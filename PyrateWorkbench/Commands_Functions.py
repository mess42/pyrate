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



from PySide.QtGui import QInputDialog
from PySide.QtGui import QLineEdit

import FreeCADGui, FreeCAD


from Object_Functions import FunctionsObject


class CreateFunctionTool:
    "Tool for creating optical system"

    def GetResources(self):
        return {"Pixmap"  : ":/icons/pyrate_func_icon.svg", # resource qrc file needed, and precompile with python-rcc
                "MenuText": "Create function ...",
                "Accel": "",
                "ToolTip": "Generates function object in document"
                }

    def IsActive(self):
        if FreeCAD.ActiveDocument == None:
            return False
        else:
            return True

    def Activated(self):

        doc = FreeCAD.ActiveDocument

        osgroups = doc.getObjectsByLabel("OS_group")
        if osgroups == []:
            FreeCAD.Console.PrintMessage("no optical system found")
        else:
            osgroup = osgroups[0]
            (name_of_functionsobject, accepted) = QInputDialog.getText(None, "Pyrate", "Name of Function Object", QLineEdit.Normal, "")        
            FunctionsObject(name_of_functionsobject, doc, osgroup) 


FreeCADGui.addCommand('CreateFunctionsCommand', CreateFunctionTool())
