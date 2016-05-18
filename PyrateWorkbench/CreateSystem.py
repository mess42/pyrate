#!/usr/bin/python3
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

        PyrateInterface.OSinterface.dummycreate4() # substitute by system creation dialog
        # dummycreate() -> lens system
        # dummycreate2() -> mirror system
        # dummycreate3() -> lens system with incorrect curvature in surface7
        # dummycreate4() -> GRIN medium
        PyrateInterface.OSinterface.createSurfaceViews(doc)
        PyrateInterface.OSinterface.showAimFiniteSurfaceStopDialog()
        PyrateInterface.OSinterface.showFieldWaveLengthDialog()
        PyrateInterface.OSinterface.createRayViews(doc, 50)
        #PyrateInterface.OSinterface.showSpotDiagrams(100)





FreeCADGui.addCommand('CreateSystemCommand', CreateSystemTool())

