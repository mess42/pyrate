#!/usr/bin/env/python
"""
Pyrate - Optical raytracing based on Python

Copyright (C) 2014-2018
               by     Moritz Esslinger moritz.esslinger@web.de
               and    Johannes Hartung j.hartung@gmx.net
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
from Interface_Identifiers import *

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
            if PyrateInterface.OSinterface.surfaceobs or PyrateInterface.OSinterface.rayobs:
                return True
            else:
                return False

    def Activated(self):

        doc = FreeCAD.ActiveDocument

        PyrateInterface.OSinterface.deleteSurfaces(doc)
        PyrateInterface.OSinterface.deleteRays(doc)

#         for so in PyrateInterface.OSinterface.surfaceobs:
#             doc.removeObject(so.Label)
#
#         for ro in PyrateInterface.OSinterface.rayobs:
#             doc.removeObject(ro.Label)
#
#         for pto in PyrateInterface.OSinterface.intersectptsobs:
#             doc.removeObject(pto.Label)
#
#         PyrateInterface.OSinterface.intersectptsobs[:] = [] # empty list
#         PyrateInterface.OSinterface.rayobs[:] = []
#         PyrateInterface.OSinterface.rayviews[:] = []
#         PyrateInterface.OSinterface.surfaceobs[:] = []
#         PyrateInterface.OSinterface.surfaceviews[:] = []

        QtGui.QMessageBox.warning(None, Title_MessageBoxes, "Notice that the optical system variable still exists and is valid!\n" + \
                                  "Only the representation of it was deleted from the active document!")

        for i in FreeCAD.ActiveDocument.Objects:
            i.touch()

        FreeCAD.ActiveDocument.recompute()


class DeleteSurfacesCommand:
    "Delete surfaces of system from active document"

    def GetResources(self):
        return {"MenuText": "Delete surfaces of system from active document ...",
                "Accel": "",
                "ToolTip": "Delete surfaces of optical system from active document",
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

        PyrateInterface.OSinterface.deleteSurfaces(doc)


        for i in FreeCAD.ActiveDocument.Objects:
            i.touch()

        FreeCAD.ActiveDocument.recompute()

class DeleteRaysCommand:
    "Delete rays through system from active document"

    def GetResources(self):
        return {"MenuText": "Delete rays through system from active document ...",
                "Accel": "",
                "ToolTip": "Delete rays through optical system from active document",
                "Pixmap": ":/icons/pyrate_del_sys_icon.svg"
                }

    def IsActive(self):
        if FreeCAD.ActiveDocument == None:
            return False
        else:
            if PyrateInterface.OSinterface.rayobs:
                return True
            else:
                return False

    def Activated(self):

        doc = FreeCAD.ActiveDocument

        PyrateInterface.OSinterface.deleteRays(doc)


        for i in FreeCAD.ActiveDocument.Objects:
            i.touch()

        FreeCAD.ActiveDocument.recompute()


FreeCADGui.addCommand('ShowRaysCommand',ShowRaysCommand())
FreeCADGui.addCommand('ShowSystemCommand',ShowSystemCommand())
FreeCADGui.addCommand('ShowSurfacesCommand',ShowSurfacesCommand())
FreeCADGui.addCommand('UpdateVisualizationCommand',UpdateVisualizationCommand())
FreeCADGui.addCommand('ShowSystemDraw2DCommand',ShowSystemDraw2DCommand())


FreeCADGui.addCommand('DeleteRaysCommand',DeleteRaysCommand())
FreeCADGui.addCommand('DeleteSystemCommand',DeleteSystemCommand())
FreeCADGui.addCommand('DeleteSurfacesCommand',DeleteSurfacesCommand())



