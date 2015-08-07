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



