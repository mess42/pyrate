import pickle
import sys
from PySide import QtGui, QtCore

import FreeCAD
import FreeCADGui
import Part
import PartGui # wichtig fuer import von icons falls keine eigenen XPMs verwendet werden
import Points


import PyrateInterface

class LoadSystemCommand:
    "Load optical system file"

    def GetResources(self):
        return {"MenuText": "Load Optical System from pickle ...",
                "Accel": "",
                "ToolTip": "Loads an optical system from pickles file",
                "Pixmap": ":/icons/pyrate_load_sys_icon.svg"
                }

    def IsActive(self):
        return True

    def Activated(self):
        if FreeCAD.ActiveDocument == None:
            FreeCAD.newDocument()


        fname, _ = QtGui.QFileDialog.getOpenFileName(None, 'Open file', '')

        if fname == "":
            return 1
        else:
            with open(fname, 'rb') as input:
                PyrateInterface.OSinterface.os = pickle.load(input)
            PyrateInterface.OSinterface.createSurfaceViews()
            #PyrateInterface.OSinterface.createRayViews()


        for i in FreeCAD.ActiveDocument.Objects:
            i.touch()

        FreeCAD.ActiveDocument.recompute()

FreeCADGui.addCommand('LoadSystemCommand',LoadSystemCommand())
