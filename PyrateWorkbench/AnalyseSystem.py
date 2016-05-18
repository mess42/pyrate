import FreeCAD
import FreeCADGui
import Part
import PartGui
import Points


import PyrateInterface
from PySide import QtGui


class ShowSpotDiagramCommand:
    "Show spot diagram"

    def GetResources(self):
        return {"MenuText": "Spot Diagram",
                "Accel": "",
                "ToolTip": "Shows the spot diagrams for all field points",
                "Pixmap": ":/icons/pyrate_del_sys_icon.svg"
                }

    def IsActive(self):
        if FreeCAD.ActiveDocument == None:
            return False
        else:
            return True

    def Activated(self):

        PyrateInterface.OSinterface.showSpotDiagrams(100)


FreeCADGui.addCommand('ShowSpotDiagramCommand',ShowSpotDiagramCommand())
