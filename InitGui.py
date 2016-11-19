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


# access to global variables for the FreeCAD interface

import FreeCADGui
import FreeCAD

from PyrateWorkbench.Interface_Checks import isLocalCoordinatesObserver

from PyrateWorkbench import Commands_OpticalSystem
from PyrateWorkbench import Observer_OpticalSystem

#from PyrateWorkbench import Commands_Files # TODO: update

from PyrateWorkbench import Commands_Visualization # TODO: update

from PyrateWorkbench import Commands_FieldPoints
from PyrateWorkbench import TaskPanel_FieldPoints

from PyrateWorkbench import Commands_Analysis # TODO: update

from PyrateWorkbench import Commands_Optimization
from PyrateWorkbench import Dialog_Optimization

from PyrateWorkbench import Commands_LocalCoordinates
from PyrateWorkbench import Observer_LocalCoordinates

from PyrateWorkbench import Commands_Functions
from PyrateWorkbench import Object_Functions


# access to the resource file
from PyrateWorkbench import resources_rc

class PyrateWorkbench ( Workbench ):
    "Pyrate workbench object"
    Icon = ":/icons/pyrate_logo_icon.svg"

    MenuText = "Pyrate Workbench"
    ToolTip = "Pyrate optical design Workbench"
    def GetClassName(self):
        return "Gui::PythonWorkbench"

    def Initialize(self):
        self.appendToolbar("Pyrate",
                           [
                           "CreateSystemCommand",
                            "Separator",
#                            "DeleteSystemCommand",
                            "UpdateVisualizationCommand",
                            "Separator",
                            "ContextAddChildToLocalCoordinatesCommand", 
                            "CreateFunctionsCommand"
#                            "LoadSystemCommand",
#                            "SaveSystemCommand"
                            ])
        #self.appendMenu("Pyrate Files", ["LoadSystemCommand", "SaveSystemCommand"]) # TODO: update
        self.appendMenu("Pyrate System", ["CreateSystemCommand", "CreateLocalCoordinatesCommand"])
        self.appendMenu("Pyrate Visualization", ["UpdateVisualizationCommand"])
        self.appendMenu("Pyrate Visualization", ["ShowSystemDraw2DCommand"])
#        self.appendMenu(["Pyrate Visualization", "Show ..."],
#                                                ["ShowSystemCommand",
#                                                 "ShowSurfacesCommand",
#                                                 "ShowRaysCommand"]
#                        )
#        self.appendMenu(["Pyrate Visualization", "Delete ..."],
#                                                ["DeleteSystemCommand",
#                                                 "DeleteSurfacesCommand",
#                                                 "DeleteRaysCommand"]
#                        )


        self.appendMenu("Pyrate Field", ["ShowAimDialogCommand", "ShowFieldDialogCommand"])
        self.appendMenu("Pyrate Analysis", ["ShowSpotDiagramCommand"])
        self.appendMenu("Pyrate Optimization", ["StartOptimizationCommand"])
        

        Log ("Loading Create System Module... done\n")

    def ContextMenu(self, recipient):
        selection = [s  for s in FreeCADGui.Selection.getSelection() if s.Document == FreeCAD.ActiveDocument ]
        Log ("selection: " + str(selection) + "\n")
        Log ("recipient: " + str(recipient) + "\n")

        if selection == []:
            self.appendContextMenu("Separator", [])
            self.appendContextMenu( "Pyrate View", 
                                       ["ShowSystemDraw2DCommand", 
                                       "ContextIncreaseScaleOfAllLocalCoordinatesCommand",
                                       "ContextDecreaseScaleOfAllLocalCoordinatesCommand"])
            self.appendContextMenu("Separator", [])
            
        
        if len(selection) == 1:
            obj = selection[0] # TODO: better classification of selections
            # TODO: why CheckObjects function not working here?
            if 'lcclass' in obj.PropertiesList:            
            #if isLocalCoordinatesObserver(obj):
                self.appendContextMenu("Separator", [])
                self.appendContextMenu( "Pyrate Local Coordinate System", 
                                       ["ContextAddChildToLocalCoordinatesCommand"])
                self.appendContextMenu("Separator", [])
            if 'wavelengths' in obj.PropertiesList:
                self.appendContextMenu("Separator", [])
                self.appendContextMenu( "Pyrate Optical System", 
                                       ["ShowFieldDialogCommand"])
                self.appendContextMenu("Separator", [])
                
                                       

    def Activated(self):
# do something here if needed...

    
        Msg ("PyrateWorkbench.Activated()\n")

    def Deactivated(self):
# do something here if needed...
        Msg ("PyrateWorkbench.Deactivated()\n")

FreeCADGui.addWorkbench(PyrateWorkbench())
