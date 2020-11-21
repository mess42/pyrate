#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Pyrate - Optical raytracing based on Python

Copyright (C) 2014-2020
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

from pprint import pformat

from pyrateoptics.raytracer.config import ConfigFile
from pyrateoptics.raytracer.material.material_glasscat import GlassCatalog

import FreeCADGui

from .Interface_Helpers import (getRelativeFilePath,
                                getAllFunctionObjects,
                                getAllMaterialCatalogues,
                                getAllLocalCoordinates,
                                getObjectByLabel)
from .Interface_Identifiers import (Material_ConstantIndexGlass,
                                    Material_ModelGlass,
                                    Material_GrinMedium,
                                    Material_CatalogMaterial1,
                                    Material_CatalogMaterial2,
                                    Material_GUI_TaskPanel_Add_TabWidget)

from .Object_Material import MaterialObject

class MaterialsTaskPanelAdd:
    def __init__(self, doc):
        # doc needs to be initialized first
        self.doc = doc

        # Catalog initialization
        databasepath = ConfigFile().get_refractive_index_database_path()
        self.glasscatalog = GlassCatalog(databasepath)

        fn = getRelativeFilePath(__file__, 'Qt/dlg_material_add.ui')
        self.form = FreeCADGui.PySideUic.loadUi(fn)

        self.form.comboBox_Shelf.clear()
        self.form.comboBox_Shelf.addItems(self.glasscatalog.get_shelves())
        self.onCurrentIndexChangedShelf(0)
        self.onCurrentIndexChangedBook(0)
        self.onCurrentIndexChangedPage(0)

        self.form.comboBox_Shelf.currentIndexChanged.connect(
            self.onCurrentIndexChangedShelf)
        self.form.comboBox_Book.currentIndexChanged.connect(
            self.onCurrentIndexChangedBook)
        self.form.comboBox_Page.currentIndexChanged.connect(
            self.onCurrentIndexChangedPage)

        self.fnobjectsindocument = getAllFunctionObjects(self.doc)

        fnobjectslabels = [obj.Label for obj in self.fnobjectsindocument]
        matcatobjectslabels = [obj.Label for obj in getAllMaterialCatalogues(doc)]
        lcobjectslabels = [obj.Label for obj in getAllLocalCoordinates(doc)]

        # this will create a Qt widget from our ui file

        self.form.comboBoxCatalogue.addItems(matcatobjectslabels)
        self.form.comboBoxLocalCoordinates.addItems(lcobjectslabels)

        # TODO: check whether glasscatalog exists and if not,
        # disable page4 and page5 of dialog

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

        self.form.pushButtonSearch.clicked.connect(self.onClickedSearchButton)
        self.form.comboBoxCandidates.currentIndexChanged.connect(
            self.onCurrentIndexChangedCandidates)

    def onClickedSearchButton(self):
        searchtext = self.form.lineEditSearch.text()
        if searchtext.strip():
            self.form.comboBoxCandidates.clear()
            self.form.comboBoxCandidates.addItems(
                [longname
                    for longname in self.glasscatalog.get_dict_of_long_names()
                    if searchtext.lower() in longname.lower()])

    def onCurrentIndexChangedCandidates(self, index):
        longname = self.form.comboBoxCandidates.currentText()
        if longname:
            mydict = {}
            shelf_book_page = self.glasscatalog.get_dict_of_long_names()[
                longname]
            mydict = self.glasscatalog.get_material_dict(*shelf_book_page)
            self.form.plainTextEditPreviewSearch.setPlainText(pformat(mydict))

    def onCurrentIndexChangedShelf(self, index):
        shelf = self.form.comboBox_Shelf.currentText()
        self.form.comboBox_Book.clear()
        if shelf != "":
            self.form.comboBox_Book.addItems(
                sorted(
                    self.glasscatalog.get_books(shelf),
                    key=lambda x: x.lower()))
            book = self.form.comboBox_Book.currentText()
            if book != "":
                self.form.comboBox_Page.clear()
                self.form.comboBox_Page.addItems(
                    sorted(self.glasscatalog.get_pages(shelf, book),
                           key=lambda x: x.lower()))

    def onCurrentIndexChangedBook(self, index):
        shelf = self.form.comboBox_Shelf.currentText()
        book = self.form.comboBox_Book.currentText()
        if shelf != "" and book != "":
            self.form.comboBox_Page.clear()
            self.form.comboBox_Page.addItems(
                sorted(self.glasscatalog.get_pages(shelf, book),
                       key=lambda x: x.lower()))

    def onCurrentIndexChangedPage(self, index):
        shelf = self.form.comboBox_Shelf.currentText()
        book = self.form.comboBox_Book.currentText()
        page = self.form.comboBox_Page.currentText()
        if shelf != "" and book != "" and page != "":
            self.form.plainTextEditPreviewChoice.setPlainText(
                pformat(self.glasscatalog.get_material_dict(shelf, book, page)))


    def updateCombo2FromCombo1WithFunctionObject(self, comb1, comb2):
        #fnobjects = self.doc.getObjectsByLabel(comb1.currentText())
        #if fnobjects != []:
        #    comb2.clear()
        #    comb2.addItems(fnobjects[0].functions)
        fnobject = getObjectByLabel(self.doc, comb1.currentText())
        if fnobject is not None:
            comb2.clear()
            comb2.addItems(fnobject.functions)


    def onCurrentIndexChangedFO_N(self, index):
        self.updateCombo2FromCombo1WithFunctionObject(
            self.form.comboBox_FO_N, self.form.comboBox_FOf_N)

    def onCurrentIndexChangedFO_Nx(self, index):
        self.updateCombo2FromCombo1WithFunctionObject(
            self.form.comboBox_FO_Nx, self.form.comboBox_FOf_Nx)

    def onCurrentIndexChangedFO_Ny(self, index):
        self.updateCombo2FromCombo1WithFunctionObject(
            self.form.comboBox_FO_Ny, self.form.comboBox_FOf_Ny)

    def onCurrentIndexChangedFO_Nz(self, index):
        self.updateCombo2FromCombo1WithFunctionObject(
            self.form.comboBox_FO_Nz, self.form.comboBox_FOf_Nz)


    # extraction functions for different material classes

    def extractConstantIndex(self):
        return {"lc": getObjectByLabel(self.doc,
                                       self.form.comboBoxLocalCoordinates.currentText()),
                "index": self.form.doubleSpinBoxIndex.value()}

    def extractModel(self):
        return {
                "lc": getObjectByLabel(self.doc,
                                       self.form.comboBoxLocalCoordinates.currentText()),
                "n0":self.form.doubleSpinBoxN0.value(),
                "A":self.form.doubleSpinBox_A.value(),
                "B":self.form.doubleSpinBox_B.value()
                }


    # TODO: GRIN: constructor changed to eat source together with function names
    def extractGrin(self):
        return {
            "fun":self.doc.getObjectsByLabel(
                self.form.comboBox_FO_N.currentText())[0].Proxy.returnSingleFunctionObject(self.form.comboBox_FOf_N.currentText()
            ),
            "dfdx":self.doc.getObjectsByLabel(
                self.form.comboBox_FO_Nx.currentText())[0].Proxy.returnSingleFunctionObject(self.form.comboBox_FOf_Nx.currentText()
            ),
            "dfdy":self.doc.getObjectsByLabel(
                self.form.comboBox_FO_Ny.currentText())[0].Proxy.returnSingleFunctionObject(self.form.comboBox_FOf_Ny.currentText()
            ),
            "dfdz":self.doc.getObjectsByLabel(
                self.form.comboBox_FO_Nz.currentText())[0].Proxy.returnSingleFunctionObject(self.form.comboBox_FOf_Nz.currentText()
            )
        }

    def extractCatalogMaterial1(self):
        return {
            "lc": getObjectByLabel(
                self.doc,
                self.form.comboBoxLocalCoordinates.currentText()),
            "ymldict": self.glasscatalog.get_material_dict(
                self.form.comboBox_Shelf.currentText(),
                self.form.comboBox_Book.currentText(),
                self.form.comboBox_Page.currentText())
            }

    def extractCatalogMaterial2(self):
        mydict = {}
        shelf_book_page = self.glasscatalog.get_dict_of_long_names()[
            self.form.comboBoxCandidates.currentText()]
        mydict = self.glasscatalog.get_material_dict(*shelf_book_page)
        return {
            "lc": getObjectByLabel(
                self.doc,
                self.form.comboBoxLocalCoordinates.currentText()),
            "ymldict": mydict
            }

    def accept(self):

        matname = self.form.lineEditName.text()

        mattype = Material_GUI_TaskPanel_Add_TabWidget[
            self.form.tabWidget.currentIndex()]

        matextractfunc = {
            Material_ConstantIndexGlass:self.extractConstantIndex,
            Material_ModelGlass:self.extractModel,
            Material_GrinMedium:self.extractGrin,
            Material_CatalogMaterial1:self.extractCatalogMaterial1,
            Material_CatalogMaterial2:self.extractCatalogMaterial2
        }

        # generating material object depending on type

        MaterialObject(self.doc,
                       self.doc.getObjectsByLabel(self.form.comboBoxCatalogue.currentText())[0],
                        matname,
                        mattype,
                        **(matextractfunc[mattype]()))

        # TODO: Check for Materials group (Label Name in Interface_Identifiers)
        # TODO: If it already exists Check for Material Mirror (if it not exists, warning, rename yes/no)
        # TODO: if there is no such group: create one with uuid name and identifier label
        # TODO: build load/save functionality


        FreeCADGui.Control.closeDialog()
