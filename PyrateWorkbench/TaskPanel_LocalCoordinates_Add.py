#!/usr/bin/env/python
"""
Pyrate - Optical raytracing based on Python

Copyright (C) 2014-2018
               by     Moritz Esslinger moritz.esslinger@web.de
               and    Johannes Hartung j.hartung@gmx.net
               and    Uwe Lippmann  uwe.lippmann@web.de
               and    Thomas Heinze t.heinze@uni-jena.de
               and    others

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

import FreeCADGui

from PySide import QtGui

from Observer_LocalCoordinates import LC

from Interface_Helpers import *
from Interface_Identifiers import *

class LocalCoordinatesTaskPanelAdd:
    def __init__(self, doc, oslabellist):
        fn = getRelativeFilePath(__file__, 'Qt/dlg_localcoords_add.ui')        
        
        # this will create a Qt widget from our ui file
        self.form = FreeCADGui.PySideUic.loadUi(fn)
        self.form.comboBoxOS.addItems(oslabellist)
        self.form.comboBoxOS.activated.connect(self.onActivatedCBOS)        

        self.doc = doc
        self.actualizeLCComboBoxFromOSComboBox()
        

    def actualizeLCComboBoxFromOSComboBox(self):
        oslabel = self.form.comboBoxOS.currentText()
        oslist = self.doc.getObjectsByLabel(oslabel)
        if oslist != []:
            osselected = oslist[0]
            lcingroup = osselected.Proxy.returnObjectsFromCoordinatesGroup()
            lclabellist = [lc.Label for lc in lcingroup]
            self.form.comboBoxParentLC.clear()
            self.form.comboBoxParentLC.addItems(lclabellist)
        else:
            self.form.comboBoxParentLC.clear()
        

    def onActivatedCBOS(self, index):
        self.actualizeLCComboBoxFromOSComboBox()        
        

    def accept(self):
        if self.form.comboBoxOS.currentText() != "":
            parentlclabel = self.form.comboBoxParentLC.currentText()        
            name_of_newlc = self.form.lineEditName.text()
            parentlclist = self.doc.getObjectsByLabel(parentlclabel)
            if parentlclist != []:
                parentlc = parentlclist[0]                
                parentlc.Proxy.addChild(name=name_of_newlc) 
            else:
                pass # either warning or new global LC
        else:
            QtGui.QMessageBox.warning(None, Title_MessageBoxes, "No optical system available! Please create one.")            


        FreeCADGui.Control.closeDialog()

        

