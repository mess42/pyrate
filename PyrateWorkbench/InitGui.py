# access to global variables for the FreeCAD interface

import FreeCADGui


import PyrateInterface
import CreateSystem
import LoadSystem
import SaveSystem
import DeleteRays
import DeleteSystem
import VisualizeSystem
import FieldOfSystem
import AnalyseSystem



# access to the resource file
import resources_rc

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
                            "DeleteSystemCommand",
                            "LoadSystemCommand",
                            "SaveSystemCommand"])
        self.appendMenu("Pyrate Files", ["LoadSystemCommand", "SaveSystemCommand"])
        self.appendMenu("Pyrate System", ["CreateSystemCommand"])
        self.appendMenu("Pyrate Visualization", ["ShowSystemCommand",
                                                 "ShowSurfacesCommand",
                                                 "ShowRaysCommand",
                                                 "DeleteSystemCommand",
                                                 "DeleteSurfacesCommand",
                                                 "DeleteRaysCommand",
                                                 "UpdateVisualizationCommand"])
        self.appendMenu("Pyrate Field", ["ShowAimDialogCommand", "ShowFieldDialogCommand"])
        self.appendMenu("Pyrate Analysis", ["ShowSpotDiagramCommand"])
        self.appendMenu("Pyrate Optimization", [""])


        Log ("Loading Create System Module... done\n")

    def Activated(self):
# do something here if needed...
        Msg ("PyrateWorkbench.Activated()\n")

    def Deactivated(self):
# do something here if needed...
        Msg ("PyrateWorkbench.Deactivated()\n")

FreeCADGui.addWorkbench(PyrateWorkbench())
