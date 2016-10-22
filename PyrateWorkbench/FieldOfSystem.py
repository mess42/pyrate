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


import time
import math
import sys
import os

import numpy as np
import matplotlib.pyplot as plt

import FreeCAD
import FreeCADGui
import Part
import PartGui
import Points


import PyrateInterface
from PySide import QtCore, QtGui


class FieldPointsTaskPanel:
    def __init__(self):
        # this will create a Qt widget from our ui file
        fn = os.path.join(os.path.dirname(__file__), 'Qt/fielddialog.ui') 
        self.form = FreeCADGui.PySideUic.loadUi(fn)
        self.form.pbAdd.clicked.connect(self.onAdd)
        self.form.pbRemove.clicked.connect(self.onRemove)
        self.form.pbLoadFile.clicked.connect(self.onLoadFile)
        self.form.pbSaveFile.clicked.connect(self.onSaveFile)

        self.form.tblFieldPoints.cellClicked.connect(self.onCellClicked)
        self.form.tblFieldPoints.cellActivated.connect(self.onCellActivated)
        
        self.row = -1        
        
    def onCellClicked(self, r, c):
        '''is cell clicked?'''
        self.row = r
    
    def onCellActivated(self, r, c):
        '''is cell activated?'''
        self.row = r

    def onAdd(self):
        '''Call Function to add field point'''
        self.form.tblFieldPoints.insertRow(self.form.tblFieldPoints.rowCount())
    
    def onRemove(self):
        '''Call Function to remove field point'''
        if self.row != -1:        
            self.form.tblFieldPoints.removeRow(self.row)
        else:
            self.form.tblFieldPoints.removeRow(self.form.tblFieldPoints.rowCount())

    def onLoadFile(self):
        '''Call Function to load field points from file'''
        pass

    def onSaveFile(self):
        '''Call Function to load field points from file'''
        pass
        


    def accept(self):
        #length = self.form.BoxLength.value()
        #width = self.form.BoxWidth.value()
        #height = self.form.BoxHeight.value()
        #if (length == 0) or (width == 0) or (height == 0):
        #    print("Error! None of the values can be 0!")
        #    # we bail out without doing anything
        #    return
        #box = Part.makeBox(length,width,height)
        #Part.show(box)
        FreeCADGui.Control.closeDialog()
        
    def reject(self):
        FreeCADGui.Control.closeDialog()
        

class ShowAimDialogCommand:
    "Show aim dialog"

    def GetResources(self):
        return {"MenuText": "Aim Dialog",
                "Accel": "",
                "ToolTip": "You may enter data for the aiming",
                "Pixmap": ":/icons/pyrate_del_sys_icon.svg"
                }

    def IsActive(self):
        if FreeCAD.ActiveDocument == None:
            return False
        else:
            return True

    def Activated(self):
        pass
        #PyrateInterface.OSinterface.showAimFiniteSurfaceStopDialog()

class ShowFieldDialogCommand:
    "Show field dialog"

    def GetResources(self):
        return {"MenuText": "Field Dialog",
                "Accel": "",
                "ToolTip": "You may enter field point data",
                "Pixmap": ":/icons/pyrate_del_sys_icon.svg"
                }

    def IsActive(self):
        if FreeCAD.ActiveDocument == None:
            return False
        else:
            return True

    def Activated(self):
        fppanel = FieldPointsTaskPanel()
        FreeCADGui.Control.showDialog(fppanel)


FreeCADGui.addCommand('ShowAimDialogCommand',ShowAimDialogCommand())
FreeCADGui.addCommand('ShowFieldDialogCommand',ShowFieldDialogCommand())



