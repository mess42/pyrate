# -*- coding: utf-8 -*-
"""
Created on Fri Nov 18 23:34:13 2016

@author: demon_ds
"""

import FreeCAD
import FreeCADGui


#from core import merit
from core import optimize

from Observer_OpticalSystem import OpticalSystemObserver
from Dialog_Optimization import OptimizationDialog

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

        # TODO: obsolete

        # non-well defined interface to internal variables of optical system
        OSinterface.os.surfaces[2].setStatus("curvature", True)
        OSinterface.os.surfaces[3].setStatus("curvature", True)
        OSinterface.os.surfaces[4].setStatus("curvature", True)
        OSinterface.os.surfaces[5].setStatus("curvature", True)
        OSinterface.os.surfaces[7].setStatus("curvature", True)

        # input

        #numsteps_t = QtGui.QLineEdit("bla", None) # 1
        #dx_t = QtGui.QLineEdit("bla", None) # 1e-6

        #FreeCAD.Console.PrintMessage(str(numsteps_t))

        optdlg = OptimizationDialog(1, 1e-6)
        optdlg.exec_()
        numsteps = optdlg.iterations
        delta = optdlg.dx

        # optimization

        #OSinterface.os = \
        #optimize.optimizeNewton1D(OSinterface.os,
        #                          merit.mySimpleDumpRMSSpotSizeMeritFunction, iterations=numsteps, dx=delta
        #                          )
        # update
        # TODO: organize in PyrateInterface class

        doc = FreeCAD.ActiveDocument

        OSinterface.deleteSurfaces(doc)
        OSinterface.deleteRays(doc)
        # abfrage!
        OSinterface.createSurfaceViews(doc)
        OSinterface.createRayViews(doc, OSinterface.shownumrays)


        for i in doc.Objects:
            i.touch()

        doc.recompute()




FreeCADGui.addCommand('StartOptimizationCommand',StartOptimizationCommand())
