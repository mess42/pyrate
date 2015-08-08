import FreeCAD
import FreeCADGui
import Part
import PartGui
import Points


import core.merit
import core.optimize

import PyrateInterface
from PySide import QtGui


class StartOptimizationCommand:
    "Starts optimization"

    def GetResources(self):
        return {"MenuText": "Start Optimization",
                "Accel": "",
                "ToolTip": "Starts Optimization",
                "Pixmap": ":/icons/pyrate_del_sys_icon.svg"
                }

    def IsActive(self):
        if FreeCAD.ActiveDocument == None:
            return False
        else:
            return True

    def Activated(self):

        # non-well defined interface to internal variables of optical system
        PyrateInterface.OSinterface.os.surfaces[2].setStatus("curvature", True)
        PyrateInterface.OSinterface.os.surfaces[3].setStatus("curvature", True)
        PyrateInterface.OSinterface.os.surfaces[4].setStatus("curvature", True)
        PyrateInterface.OSinterface.os.surfaces[5].setStatus("curvature", True)
        PyrateInterface.OSinterface.os.surfaces[7].setStatus("curvature", True)

        # optimization

        PyrateInterface.OSinterface.os = \
        core.optimize.optimizeNewton1D(
                                       PyrateInterface.OSinterface.os,
                                       core.merit.myPersonalMeritFunctionForTestingPurposes, iterations=1, dx=1e-6
                                       )
        # update (todo: organize in PyrateInterface class)

        doc = FreeCAD.ActiveDocument

        PyrateInterface.OSinterface.deleteSurfaces(doc)
        PyrateInterface.OSinterface.deleteRays(doc)
        # abfrage!
        PyrateInterface.OSinterface.createSurfaceViews(doc)
        PyrateInterface.OSinterface.createRayViews(doc, PyrateInterface.OSinterface.shownumrays)


        for i in doc.Objects:
            i.touch()

        doc.recompute()




FreeCADGui.addCommand('StartOptimizationCommand',StartOptimizationCommand())
