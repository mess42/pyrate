#!/usr/bin/python3
import pickle
from PySide import QtGui, QtCore

import FreeCAD
import FreeCADGui
import Part
import PartGui # wichtig fuer import von icons falls keine eigenen XPMs verwendet werden
import Points


import PyrateInterface

class SaveSystemCommand:
    "Save optical system file"

    def GetResources(self):
        return {"MenuText": "Save Optical System as pickle ...",
                "Accel": "",
                "ToolTip": "Saves an optical system to pickles file",
                "Pixmap": ":/icons/pyrate_save_sys_icon.svg"
                }

    def IsActive(self):
        if FreeCAD.ActiveDocument == None:
            return False
        else:
            return True

    def Activated(self):
        savedlg = QtGui.QFileDialog(None, 'Save file', '')
        savedlg.setFileMode(QtGui.QFileDialog.AnyFile)
        fname, _ = savedlg.getSaveFileName()


        if fname == "":
            return 1
        else:
            with open(fname, 'wb') as output:
                pickle.dump(PyrateInterface.OSinterface.os, output, pickle.HIGHEST_PROTOCOL)
            QtGui.QMessageBox.warning(None, "Pyrate", "Dumped!\n" + "Warning, pickle format may change over time!")


FreeCADGui.addCommand('SaveSystemCommand',SaveSystemCommand())
