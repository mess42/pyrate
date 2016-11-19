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
import os
import FreeCADGui

from PySide.QtGui import QInputDialog, QLineEdit

from Object_Functions import FunctionsObject

class FunctionsTaskPanelAdd:
    def __init__(self, doc, stringlist):
        fn = os.path.join(os.path.dirname(__file__), 'Qt/dlg_functionsobject_add.ui')        
        
        # this will create a Qt widget from our ui file
        self.form = FreeCADGui.PySideUic.loadUi(fn)
        self.form.comboBox.addItems(stringlist)
        self.doc = doc

    def accept(self):
        oslabel = self.form.comboBox.currentText()
        name_of_functionsobject = self.form.lineEditName.text()
        
        os = self.doc.getObjectsByLabel(oslabel)[0]
                
        fngroupname = os.NameFunctionsGroup
        fngroup = self.doc.getObject(fngroupname)

        FunctionsObject(name_of_functionsobject, self.doc, fngroup) 

        FreeCADGui.Control.closeDialog()

        

