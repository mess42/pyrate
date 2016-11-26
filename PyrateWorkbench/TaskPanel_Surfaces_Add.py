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
import FreeCADGui

from Interface_Helpers import *
from Interface_Checks import *

class SurfacesTaskPanelAdd:
    def __init__(self, doc):
        # doc needs to be initialized first        
        self.doc = doc


        fn = getRelativeFilePath(__file__, 'Qt/dlg_surface_add.ui')        
        self.form = FreeCADGui.PySideUic.loadUi(fn)

        # now add all optical systems to combobox
        self.form.comboBoxOS.activated.connect(self.onActivatedCBOS)        


        self.form.comboBoxOS.clear()
        self.form.comboBoxOS.addItems([os.Label for os in getAllOpticalSystemObservers(self.doc)])
        
        self.form.comboBoxLC.clear()
        self.updateComboLCfromComboOS(self.form.comboBoxOS, self.form.comboBoxLC)
        
        self.form.comboBoxMaterial.clear()
        self.form.comboBoxMaterial.addItems([mat.Label for mat in getAllMaterials(self.doc)])

    def updateComboLCfromComboOS(self, comboos, combolc):
        labelos = comboos.currentText()
        oss = self.doc.getObjectsByLabel(labelos)
        if oss != None:
            combolc.clear()
            combolc.addItems([lc.Label for lc in oss[0].Proxy.returnObjectsFromCoordinatesGroup()])
    
    def onActivatedCBOS(self, index):
        self.updateComboLCfromComboOS(self.form.comboBoxOS, self.form.comboBoxLC)

    def accept(self):
        

        FreeCADGui.Control.closeDialog()
