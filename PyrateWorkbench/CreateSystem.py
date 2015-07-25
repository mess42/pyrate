import FreeCAD,FreeCADGui, Part
import PartGui # wichtig fuer import von icons falls keine eigenen XPMs verwendet werden
import Points

import pickle

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


class LoadSystemCommand:
    "Load optical system file"

    def GetResources(self):
        return {"MenuText": "Load Optical System ...",
                "Accel": "Ctrl+L",
                "ToolTip": "Loads an optical system from pickles file",
                "Pixmap": ":/icons/File_Document-open.svg"
#                "Pixmap": ":/icons/Part_Point_Parametric.svg" # standard icon aus FreeCAD
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
        return True

    def Activated(self):
#        PyrateInterface.OSinterface.configure_traits()
#        FreeCAD.Console.PrintMessage(str(PyrateInterface.OSinterface.bla) + "\n")

        if FreeCAD.ActiveDocument == None:
            # create new document
            return False


        fname, _ = QtGui.QFileDialog.getOpenFileName(None, 'Open file', '/home')


        if fname == "":
            return 1
        else:
            with open(fname, 'rb') as input:
                PyrateInterface.OSinterface.os = pickle.load(input)
            PyrateInterface.OSinterface.createSurfaceViews()
            PyrateInterface.OSinterface.createRayViews()


        for i in FreeCAD.ActiveDocument.Objects:
            i.touch()

        FreeCAD.ActiveDocument.recompute()





FreeCADGui.addCommand('CreateSystemCommand', CreateSystemTool())
FreeCADGui.addCommand('LoadSystemCommand',LoadSystemCommand())
