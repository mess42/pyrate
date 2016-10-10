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

from PyrateWorkbench import PyrateInterface
from PyrateWorkbench import CreateSystem
from PyrateWorkbench import LoadSystem
from PyrateWorkbench import SaveSystem
from PyrateWorkbench import DeleteRays
from PyrateWorkbench import DeleteSystem
from PyrateWorkbench import VisualizeSystem
from PyrateWorkbench import FieldOfSystem
from PyrateWorkbench import AnalyseSystem
from PyrateWorkbench import OptimizeSystem

from PyrateWorkbench import LocalCoordinatesTree

#import PyrateWorkbench.PyrateInterface
#import PyrateWorkbench.CreateSystem
#import PyrateWorkbench.LoadSystem
#import PyrateWorkbench.SaveSystem
#import PyrateWorkbench.DeleteRays
#import PyrateWorkbench.DeleteSystem
#import PyrateWorkbench.VisualizeSystem
#import PyrateWorkbench.FieldOfSystem
#import PyrateWorkbench.AnalyseSystem
#import PyrateWorkbench.OptimizeSystem



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
                           ["CreateSystemCommand",
                            "Separator",
                            "DeleteSystemCommand",
                            "UpdateVisualizationCommand",
                            "Separator",
                            "LoadSystemCommand",
                            "SaveSystemCommand"])
        self.appendMenu("Pyrate Files", ["LoadSystemCommand", "SaveSystemCommand"])
        self.appendMenu("Pyrate System", ["CreateSystemCommand", "CreateLocalCoordinatesCommand"])
        self.appendMenu("Pyrate Visualization", ["UpdateVisualizationCommand"])
        self.appendMenu("Pyrate Visualization", ["ShowSystemDraw2DCommand"])
        self.appendMenu(["Pyrate Visualization", "Show ..."],
                                                ["ShowSystemCommand",
                                                 "ShowSurfacesCommand",
                                                 "ShowRaysCommand"]
                        )
        self.appendMenu(["Pyrate Visualization", "Delete ..."],
                                                ["DeleteSystemCommand",
                                                 "DeleteSurfacesCommand",
                                                 "DeleteRaysCommand"]
                        )


        self.appendMenu("Pyrate Field", ["ShowAimDialogCommand", "ShowFieldDialogCommand"])
        self.appendMenu("Pyrate Analysis", ["ShowSpotDiagramCommand"])
        self.appendMenu("Pyrate Optimization", ["StartOptimizationCommand"])


        Log ("Loading Create System Module... done\n")

    def Activated(self):
# do something here if needed...
        Msg ("PyrateWorkbench.Activated()\n")

    def Deactivated(self):
# do something here if needed...
        Msg ("PyrateWorkbench.Deactivated()\n")

FreeCADGui.addWorkbench(PyrateWorkbench())
