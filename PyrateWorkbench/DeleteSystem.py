import FreeCAD
import FreeCADGui
import Part
import PartGui # wichtig fuer import von icons falls keine eigenen XPMs verwendet werden
import Points


import PyrateInterface
from PySide import QtGui

class DeleteSystemCommand:
    "Delete system from active document"

    def GetResources(self):
        return {"MenuText": "Delete optical system from active document ...",
                "Accel": "",
                "ToolTip": "Delete optical system from active document (including rays and intersection points)",
                "Pixmap": ":/icons/pyrate_del_sys_icon.svg"
                }

    def IsActive(self):
        if FreeCAD.ActiveDocument == None:
            return False
        else:
            if PyrateInterface.OSinterface.surfaceobs:
                return True
            else:
                return False

    def Activated(self):

        doc = FreeCAD.ActiveDocument

        for so in PyrateInterface.OSinterface.surfaceobs:
            doc.removeObject(so.Label)

        for ro in PyrateInterface.OSinterface.rayobs:
            doc.removeObject(ro.Label)

        for pto in PyrateInterface.OSinterface.intersectptsobs:
            doc.removeObject(pto.Label)

        PyrateInterface.OSinterface.intersectptsobs[:] = [] # empty list
        PyrateInterface.OSinterface.rayobs[:] = []
        PyrateInterface.OSinterface.rayviews[:] = []
        PyrateInterface.OSinterface.surfaceobs[:] = []
        PyrateInterface.OSinterface.surfaceviews[:] = []

        QtGui.QMessageBox.warning(None, "Pyrate", "Notice that the optical system variable still exists and is valid!\n" + \
                                  "Only the representation of it was deleted from the active document!")

        for i in FreeCAD.ActiveDocument.Objects:
            i.touch()

        FreeCAD.ActiveDocument.recompute()

FreeCADGui.addCommand('DeleteSystemCommand',DeleteSystemCommand())
