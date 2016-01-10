#!/usr/bin/env/python
"""
Pyrate - Optical raytracing based on Python

Copyright (C) 2014 Moritz Esslinger moritz.esslinger@web.de
               and Johannes Hartung j.hartung@gmx.net
               and    Uwe Lippmann  uwe.lippmann@web.de

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 2
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
"""


import FreeCAD
import FreeCADGui
import Part
import PartGui
import Points


from core import merit
from core import optimize

import PyrateInterface
from PySide import QtGui


class OptimizeDialog(QtGui.QDialog):
    def __init__(self, iters, delta):
        super(OptimizeDialog, self).__init__()

        self.iterations = iters
        self.dx = delta

        self.initUI()

    def initUI(self):

        lbliters = QtGui.QLabel('Iterations', self)
        lbliters.move(10,10)

        self.qliters = QtGui.QLineEdit(str(self.iterations), self)
        self.qliters.move(10, 50)

        lbliters = QtGui.QLabel('dx', self)
        lbliters.move(10,90)

        self.qldx = QtGui.QLineEdit(str(self.dx), self)
        self.qldx.move(10, 130)

        okbtn = QtGui.QPushButton("OK", self)
        okbtn.move(10, 170)
        okbtn.clicked.connect(self.onOK)


        self.setGeometry(300, 300, 600, 500)
        self.setWindowTitle('Field Configuration Dialog')
        self.show()

    def onOK(self):
        self.iters = int(self.qliters.text())
        self.dx = float(self.qldx.text())
        self.close()


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

        # input

        #numsteps_t = QtGui.QLineEdit("bla", None) # 1
        #dx_t = QtGui.QLineEdit("bla", None) # 1e-6

        #FreeCAD.Console.PrintMessage(str(numsteps_t))

        optdlg = OptimizeDialog(1, 1e-6)
        optdlg.exec_()
        numsteps = optdlg.iterations
        delta = optdlg.dx

        # optimization

        PyrateInterface.OSinterface.os = \
        optimize.optimizeNewton1D(PyrateInterface.OSinterface.os,
                                  merit.mySimpleDumpRMSSpotSizeMeritFunction, iterations=numsteps, dx=delta
                                  )
        # update
        # TODO: organize in PyrateInterface class

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
