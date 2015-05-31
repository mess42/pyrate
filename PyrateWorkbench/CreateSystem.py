import FreeCAD,FreeCADGui, Part
import PartGui # wichtig fuer import von icons falls keine eigenen XPMs verwendet werden
import Points

from traits.api import HasTraits, Str, Int, Float
import traitsui

from PySide import QtGui, QtCore

from numpy import *

import os.path
ICONS_PATH = os.path.dirname(__file__)

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
        QtGui.QMessageBox.information(None,"","Houston, we have a problem")
        ball = Ball()
        ball.configure_traits()

        pp = Points.Points()
        ptslist =[]
        for i in range(1000):
            r3 = random.randn(3)
            r3 = ball.radius * r3/linalg.norm(r3)
            ptslist.append(tuple(r3))
        pp.addPoints(ptslist)
        Points.show(pp)
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
        r1 = random.rand(3) + 0.1
        r2 = 2.0*random.rand(3) - 1.0
        r3 = random.randn(3)
        r3 = r3/linalg.norm(r3)
        FreeCAD.Console.PrintMessage("rand1: " + str(r1) + " rand2: " + str(r2) + "\n")
        maxl = 10.0;
        maxrange = 20;
        cube = Part.makeBox(maxl*r1[0], maxl*r1[1], maxl*r1[2], FreeCAD.Vector(tuple(maxrange*r2)), FreeCAD.Vector(tuple(r3)))
        Part.show(cube)

        #f=ray.Wrapper(ray.Ray,[0,0,0],[1,2,3])
        #FreeCAD.Console.PrintMessage(f.__doc__)
        #FreeCAD.Console.PrintMessage(str(f.getStart()))
        #FreeCAD.Console.PrintMessage(str(f.getEnd()))

        for i in FreeCAD.ActiveDocument.Objects:
            i.touch()

        FreeCAD.ActiveDocument.recompute()





FreeCADGui.addCommand('CreateSystemCommand', CreateSystemTool())
#FreeCADGui.addCommand('MyCommand2',MyTool2())
