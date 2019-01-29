#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 29 23:03:36 2019

@author: joha2
"""

import FreeCAD, FreeCADGui
from TaskPanel_ClassWithOptimizableVariables_Edit import ClassWithOptimizableVariablesTaskPanelEdit

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

