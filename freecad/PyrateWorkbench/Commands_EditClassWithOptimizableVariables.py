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

import FreeCAD, FreeCADGui
from .TaskPanel_ClassWithOptimizableVariables_Edit import ClassWithOptimizableVariablesTaskPanelEdit

class EditClassWithOptimizableVariablesTool:
    "Tool for editing class with optimizable variables object"

    def GetResources(self):
        return {"Pixmap"  : ":/icons/pyrate_material_catalogue_icon.svg", # resource qrc file needed, and precompile with python-rcc
                "MenuText": "Edit class with optimizable variables",
                "Accel": "",
                "ToolTip": "Edit class with optimizable variable in document"
                }

    def IsActive(self):
        return True

    def Activated(self):
        cwov_panel = ClassWithOptimizableVariablesTaskPanelEdit(None)
        # wrong call (can be set to an object_reference from tree)
        FreeCADGui.Control.showDialog(cwov_panel)

FreeCADGui.addCommand('EditClassWithOptimizableVariablesCommand',
                      EditClassWithOptimizableVariablesTool())

