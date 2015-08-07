from PySide import QtGui, QtCore

import FreeCAD
import FreeCADGui
import Part
import PartGui # wichtig fuer import von icons falls keine eigenen XPMs verwendet werden
import Points


import PyrateInterface


class CreateSystemTool:
    "Tool for creating optical system"

    def GetResources(self):
        return {"Pixmap"  : ":/icons/pyrate_logo_icon.svg", # resource qrc file needed, and precompile with python-rcc
                "MenuText": "Create optical system ...",
                "Accel": "",
                "ToolTip": "Opens dialog for system creation"
                }

    def IsActive(self):
        if FreeCAD.ActiveDocument == None:
            return False
        else:
            return True

    def Activated(self):

        doc = FreeCAD.ActiveDocument

        PyrateInterface.OSinterface.dummycreate() # substitute by system creation dialog
        PyrateInterface.OSinterface.createSurfaceViews(doc)
        PyrateInterface.OSinterface.showAimFiniteSurfaceStopDialog()
        PyrateInterface.OSinterface.showFieldWaveLengthDialog()
        PyrateInterface.OSinterface.createRayViews(doc, 10)
        PyrateInterface.OSinterface.showSpotDiagrams(100)





FreeCADGui.addCommand('CreateSystemCommand', CreateSystemTool())

