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

class MaterialsTaskPanelAdd:
    def __init__(self, doc):
        # doc needs to be initialized first        
        self.doc = doc


        fn = getRelativeFilePath(__file__, 'Qt/dlg_material_add.ui')        

        fnobjectsindocument = []
        for obj in doc.Objects:
           if isFunctionsObject(obj):
               fnobjectsindocument.append(obj)
               
        fnobjectslabels = [obj.Label for obj in fnobjectsindocument]
        matcatobjectslabels = [obj.Label for obj in getAllMaterialCatalogues(doc)]
        
        # this will create a Qt widget from our ui file
        self.form = FreeCADGui.PySideUic.loadUi(fn)
        
        self.form.comboBoxCatalogue.addItems(matcatobjectslabels)
                
        
        
        self.form.comboBox_FO_N.addItems(fnobjectslabels)
        self.form.comboBox_FO_Nx.addItems(fnobjectslabels)
        self.form.comboBox_FO_Ny.addItems(fnobjectslabels)
        self.form.comboBox_FO_Nz.addItems(fnobjectslabels)
        
        self.form.comboBox_FO_N.currentIndexChanged.connect(self.onCurrentIndexChangedFO_N)
        self.form.comboBox_FO_Nx.currentIndexChanged.connect(self.onCurrentIndexChangedFO_Nx)        
        self.form.comboBox_FO_Ny.currentIndexChanged.connect(self.onCurrentIndexChangedFO_Ny)        
        self.form.comboBox_FO_Nz.currentIndexChanged.connect(self.onCurrentIndexChangedFO_Nz)        

        self.updateCombo2FromCombo1WithFunctionObject(self.form.comboBox_FO_N, self.form.comboBox_FOf_N)
        self.updateCombo2FromCombo1WithFunctionObject(self.form.comboBox_FO_Nx, self.form.comboBox_FOf_Nx)
        self.updateCombo2FromCombo1WithFunctionObject(self.form.comboBox_FO_Ny, self.form.comboBox_FOf_Ny)
        self.updateCombo2FromCombo1WithFunctionObject(self.form.comboBox_FO_Nz, self.form.comboBox_FOf_Nz)
        
        

    def updateCombo2FromCombo1WithFunctionObject(self, comb1, comb2):
        fnobjects = self.doc.getObjectsByLabel(comb1.currentText())
        if fnobjects != []:
            comb2.clear()
            comb2.addItems(fnobjects[0].functions)
        

    def onCurrentIndexChangedFO_N(self, index):
        self.updateCombo2FromCombo1WithFunctionObject(self.form.comboBox_FO_N, self.form.comboBox_FOf_N)

    def onCurrentIndexChangedFO_Nx(self, index):
        self.updateCombo2FromCombo1WithFunctionObject(self.form.comboBox_FO_Nx, self.form.comboBox_FOf_Nx)

    def onCurrentIndexChangedFO_Ny(self, index):
        self.updateCombo2FromCombo1WithFunctionObject(self.form.comboBox_FO_Ny, self.form.comboBox_FOf_Ny)

    def onCurrentIndexChangedFO_Nz(self, index):
        self.updateCombo2FromCombo1WithFunctionObject(self.form.comboBox_FO_Nz, self.form.comboBox_FOf_Nz)


    def accept(self):
        # TODO: Check for Materials group (Label Name in Interface_Identifiers)
        # TODO: If it already exists Check for Material Mirror (if it not exists, warning, rename yes/no)        
        # TODO: if there is no such group: create one with uuid name and identifier label
        # TODO: build load/save functionality
        
        
        #oslabel = self.form.comboBox.currentText()
        #name_of_functionsobject = self.form.lineEditName.text()
        
        #os = self.doc.getObjectsByLabel(oslabel)[0]
        #        
        #fngroupname = os.NameFunctionsGroup
        #fngroup = self.doc.getObject(fngroupname)

        #FunctionsObject(name_of_functionsobject, self.doc, fngroup) 

        

        FreeCADGui.Control.closeDialog()
