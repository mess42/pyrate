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

from Object_Functions import FunctionsObject, FunctionsView

from Interface_Helpers import *
from Interface_Identifiers import *

class FunctionsTaskPanelAdd:
    def __init__(self, doc, stringlist):
        fn = getRelativeFilePath(__file__, 'Qt/dlg_functionsobject_add.ui')        
        
        # this will create a Qt widget from our ui file
        self.form = FreeCADGui.PySideUic.loadUi(fn)
        self.form.comboBox.addItems(stringlist)
        self.doc = doc

    def accept(self):
        oslabel = self.form.comboBox.currentText()
        name_of_functionsobject = self.form.lineEditName.text()
        initial_source_code = self.form.plainTextEdit.toPlainText()
        
        try:        
            os = self.doc.getObjectsByLabel(oslabel)[0]
        except IndexError:
            QtGui.QMessageBox.warning(None, Title_MessageBoxes, "No optical system available! Please create one.")            
        else:
            fngroupname = os.NameFunctionsGroup
            fngroup = self.doc.getObject(fngroupname)
    
            fnobj = FunctionsObject(name_of_functionsobject, initial_source_code, self.doc, fngroup) 
            FunctionsView(fnobj.Object.ViewObject)

        FreeCADGui.Control.closeDialog()

        

