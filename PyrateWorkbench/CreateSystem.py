import FreeCAD,FreeCADGui, Part
import PartGui # wichtig fuer import von icons falls keine eigenen XPMs verwendet werden
import Points

from traits.api import HasTraits, Str, Int, Float
import traitsui

from PySide import QtGui, QtCore

from numpy import *

import PyrateInterface


class Ball(HasTraits):
    radius = Float

class CreateSystemTool:
    "Tool for creating optical system"

    def GetResources(self):
        return {"Pixmap"  : ":/icons/pyrate_logo_icon.svg", # resource qrc file needed, and precompile with python-rcc
                "MenuText": "Create optical system ...",
                "Accel": "Ctrl+M",
                "ToolTip": "Opens dialog for system creation"
                }

    def IsActive(self):
        if FreeCAD.ActiveDocument == None:
            return False
        else:
            return True

    def Activated(self):
        PyrateInterface.OSinterface.dummycreate2()
        PyrateInterface.OSinterface.createSurfaceViews()
        PyrateInterface.OSinterface.createRayViews()
        # would be useful if one could center the view on a specific surface
        # how to implement apertures?
        # todo: create view objects for rays and for intersection points
        # draw coordinate frame per default
        # selection coupling to certain surface
        # later create system by table and update view per button

#        QtGui.QMessageBox.information(None,"","Houston, we have a problem")
#        ball = Ball()
#        ball.configure_traits()

#        pp = Points.Points()
#        ptslist =[]
#        for i in range(1000):
#            r3 = random.randn(3)
#            r3 = ball.radius * r3/linalg.norm(r3)
#            ptslist.append(tuple(r3))
#        pp.addPoints(ptslist)
#        Points.show(pp)
# do something here...


class MyTool2:
    "My tool2 object"

    def GetResources(self):
        return {"MenuText": "My Command 2",
                "Accel": "Ctrl+N",
                "ToolTip": "My extraordinary command2",
                "Pixmap": ":/icons/Part_Point_Parametric.svg" # standard icon aus FreeCAD
##                """
##                /* XPM */
##                static const char *test_icon2[]={
##                "16 16 2 1",
##                "a c #000000",
##                ". c None",
##                "................",
##                "................",
##                "..############..",
##                "..############..",
##                "..############..",
##                "......####......",
##                "......####......",
##                "......####......",
##                "......####......",
##                "......####......",
##                "......####......",
##                "......####......",
##                "......####......",
##                "......####......",
##                "................",
##                "................"};
##                """
                }

    def IsActive(self):
        if FreeCAD.ActiveDocument == None:
            return False
        else:
            return True

    def Activated(self):
        PyrateInterface.OSinterface.configure_traits()
        FreeCAD.Console.PrintMessage(str(PyrateInterface.OSinterface.bla) + "\n")

        for i in FreeCAD.ActiveDocument.Objects:
            i.touch()

        FreeCAD.ActiveDocument.recompute()





FreeCADGui.addCommand('CreateSystemCommand', CreateSystemTool())
FreeCADGui.addCommand('MyCommand2',MyTool2())
