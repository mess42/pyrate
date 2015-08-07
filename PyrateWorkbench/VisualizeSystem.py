import FreeCAD
import FreeCADGui
import Part
import PartGui
import Points


import PyrateInterface
from PySide import QtGui


class ShowSystemCommand:
    "Show system in active document"

    def GetResources(self):
        return {"MenuText": "Show complete system in active document",
                "Accel": "",
                "ToolTip": "Show all system components in active document",
                "Pixmap": ":/icons/pyrate_del_sys_icon.svg"
                }

    def IsActive(self):
        if FreeCAD.ActiveDocument == None:
            return False
        else:
            if not PyrateInterface.OSinterface.surfaceobs and not PyrateInterface.OSinterface.rayobs:
                return True
            else:
                return False

    def Activated(self):

        doc = FreeCAD.ActiveDocument

        PyrateInterface.OSinterface.createSurfaceViews(doc)
        # ask for rayviews
        PyrateInterface.OSinterface.createRayViews(10, doc)

        for i in doc.Objects:
            i.touch()

        doc.recompute()


class ShowSurfacesCommand:
    "Show surfaces of system in active document"

    def GetResources(self):
        return {"MenuText": "Show surfaces in active document",
                "Accel": "",
                "ToolTip": "Show surfaces in active document",
                "Pixmap": ":/icons/pyrate_del_sys_icon.svg"
                }

    def IsActive(self):
        if FreeCAD.ActiveDocument == None:
            return False
        else:
            if not PyrateInterface.OSinterface.surfaceobs:
                return True
            else:
                return False

    def Activated(self):

        doc = FreeCAD.ActiveDocument

        PyrateInterface.OSinterface.createSurfaceViews(doc)


        for i in doc.Objects:
            i.touch()

        doc.recompute()

class ShowRaysCommand:
    "Delete surfaces of system from active document"

    def GetResources(self):
        return {"MenuText": "Show rays through system in active document",
                "Accel": "",
                "ToolTip": "Show rays through optical system in active document",
                "Pixmap": ":/icons/pyrate_del_sys_icon.svg"
                }

    def IsActive(self):
        if FreeCAD.ActiveDocument == None:
            return False
        else:
            if not PyrateInterface.OSinterface.rayobs:
                return True
            else:
                return False

    def Activated(self):

        doc = FreeCAD.ActiveDocument

        # abfrage!
        PyrateInterface.OSinterface.createRayViews(doc, 10)


        for i in doc.Objects:
            i.touch()

        doc.recompute()


class UpdateVisualizationCommand:
    "Update System representation in active document"

    def GetResources(self):
        return {"MenuText": "Update System in active document",
                "Accel": "Ctrl+U",
                "ToolTip": "Updates representation of optical system in active document",
                "Pixmap": ":/icons/pyrate_del_sys_icon.svg"
                }

    def IsActive(self):
        if FreeCAD.ActiveDocument == None:
            return False
        else:
            return True

    def Activated(self):

        doc = FreeCAD.ActiveDocument

        PyrateInterface.OSinterface.deleteSurfaces(doc)
        PyrateInterface.OSinterface.deleteRays(doc)
        # abfrage!
        PyrateInterface.OSinterface.createSurfaceViews(doc)
        PyrateInterface.OSinterface.createRayViews(doc, 10)


        for i in doc.Objects:
            i.touch()

        doc.recompute()



FreeCADGui.addCommand('ShowRaysCommand',ShowRaysCommand())
FreeCADGui.addCommand('ShowSystemCommand',ShowSystemCommand())
FreeCADGui.addCommand('ShowSurfacesCommand',ShowSurfacesCommand())
FreeCADGui.addCommand('UpdateVisualizationCommand',UpdateVisualizationCommand())

