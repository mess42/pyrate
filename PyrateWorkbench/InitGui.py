# access to global variables for the FreeCAD interface

import FreeCADGui


import PyrateInterface
import CreateSystem
import LoadSystem
import SaveSystem
import DeleteRays
import DeleteSystem



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
        self.appendToolbar("Pyrate", ["CreateSystemCommand", "LoadSystemCommand", "SaveSystemCommand"])
        self.appendMenu("Pyrate Main", ["CreateSystemCommand"])
        self.appendMenu("Pyrate Files", ["LoadSystemCommand", "SaveSystemCommand"])

        Log ("Loading Create System Module... done\n")

    def Activated(self):
# do something here if needed...
        Msg ("PyrateWorkbench.Activated()\n")

    def Deactivated(self):
# do something here if needed...
        Msg ("PyrateWorkbench.Deactivated()\n")

FreeCADGui.addWorkbench(PyrateWorkbench())
