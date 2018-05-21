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

import platform
import os
import sys

# TODO: at the moment only available for scipy, add further libs in the future
try:
    print("checking for scipy support")
    import scipy
except:
    print("scipy not found:")
    print("please add appropriate paths to paths.txt, i.e. c:/yourpython/Lib/site-packages/")
    if platform.system() == 'Windows':

        winpathfile = os.path.join(os.path.dirname(__file__), "paths.txt")
        wpfile = open(winpathfile)
        for line in wpfile:
            sys.path.append(line)
        wpfile.close()
    if platform.system() == 'Linux':
        print("please install scipy via your favourite package manager or download an appropriate python distribution")
    print("retrying to import")
    import scipy
else:
    print("found scipy")
    


# access to global variables for the FreeCAD interface

import FreeCADGui
import FreeCAD

from freecad.PyrateWorkbench.Interface_Checks import isLocalCoordinatesObserver

from freecad.PyrateWorkbench import Commands_OpticalSystem
from freecad.PyrateWorkbench import Observer_OpticalSystem

#from freecad.PyrateWorkbench import Commands_Files # TODO: update

from freecad.PyrateWorkbench import Commands_Visualization # TODO: update

from freecad.PyrateWorkbench import Commands_FieldPoints
from freecad.PyrateWorkbench import TaskPanel_FieldPoints

from freecad.PyrateWorkbench import Commands_Analysis # TODO: update

from freecad.PyrateWorkbench import Commands_Optimization
from freecad.PyrateWorkbench import Dialog_Optimization

from freecad.PyrateWorkbench import Commands_LocalCoordinates
from freecad.PyrateWorkbench import Observer_LocalCoordinates

from freecad.PyrateWorkbench import Commands_Functions
from freecad.PyrateWorkbench import Object_Functions

from freecad.PyrateWorkbench import Commands_Materials

from freecad.PyrateWorkbench import Commands_Surface

# access to the resource file
from freecad.PyrateWorkbench import resources_rc

class PyrateWorkbench ( FreeCADGui.Workbench ):
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
#                            "UpdateVisualizationCommand",
#                            "Separator",
                            "CreateLocalCoordinatesCommand", 
                            "CreateFunctionsCommand",
                            "Separator",
                            "CreateMaterialsCatalogueCommand",
                            "CreateMaterialsCommand",
                            "Separator",
                            "CreateSurfacesCommand"
#                            "LoadSystemCommand",
#                            "SaveSystemCommand"
                            ])
        #self.appendMenu("Pyrate Files", ["LoadSystemCommand", "SaveSystemCommand"]) # TODO: update
        self.appendMenu("Pyrate System", 
                        ["CreateSystemCommand", 
                        "CreateLocalCoordinatesCommand",
                        "CreateFunctionsCommand",
                        "CreateSurfacesCommand",
                        "ShowSurfaceDialogCommand"
                        ])
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


        self.appendMenu("Pyrate Field", 
                        [
#                        "ShowAimDialogCommand", 
                        "ShowFieldDialogCommand"
                        ])
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
                self.appendContextMenu( "Pyrate Field", 
                                       ["ShowFieldDialogCommand"])
                self.appendContextMenu( "Pyrate Surfaces", 
                                       ["ShowSurfaceDialogCommand"])
                self.appendContextMenu("ShowRaybundlesCommand", [])
                self.appendContextMenu("Separator", [])
                
                                       

    def Activated(self):
# do something here if needed...

    
        Msg ("PyrateWorkbench.Activated()\n")

    def Deactivated(self):
# do something here if needed...
        Msg ("PyrateWorkbench.Deactivated()\n")

FreeCADGui.addWorkbench(PyrateWorkbench())
