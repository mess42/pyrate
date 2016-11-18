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


from Observer_OpticalSystem import OpticalSystemObserver
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
        PyrateInterface.OSinterface.createRayViews(PyrateInterface.OSinterface.shownumrays, doc)

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
        PyrateInterface.OSinterface.createRayViews(doc, PyrateInterface.OSinterface.shownumrays)


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
        PyrateInterface.OSinterface.createRayViews(doc, PyrateInterface.OSinterface.shownumrays)


        for i in doc.Objects:
            i.touch()

        doc.recompute()


class ShowSystemDraw2DCommand:
    "Show optical system in draw2d perspective"

    def GetResources(self):
        return {"MenuText": "Draw2d perspective",
                "Accel": "",
                "ToolTip": "Shows optical system in draw2d representation",
                "Pixmap": ":/icons/pyrate_del_sys_icon.svg"
                }

    def IsActive(self):
        if FreeCAD.ActiveDocument == None:
            return False
        else:
            return True

    def Activated(self):

        docview = FreeCADGui.ActiveDocument.ActiveView

        docview.viewLeft()
        docview.viewRotateRight()



FreeCADGui.addCommand('ShowRaysCommand',ShowRaysCommand())
FreeCADGui.addCommand('ShowSystemCommand',ShowSystemCommand())
FreeCADGui.addCommand('ShowSurfacesCommand',ShowSurfacesCommand())
FreeCADGui.addCommand('UpdateVisualizationCommand',UpdateVisualizationCommand())
FreeCADGui.addCommand('ShowSystemDraw2DCommand',ShowSystemDraw2DCommand())

