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
