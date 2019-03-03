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

from pprint import pprint

import FreeCADGui
from PySide.QtGui import QTableWidgetItem, QHeaderView
from PySide.QtCore import Qt

from pyrateoptics.core.base_ui import UIInterfaceClassWithOptimizableVariables

from .Interface_Helpers import getRelativeFilePath


class ClassWithOptimizableVariablesTaskPanelEdit:
    """
    User interface to edit class with optimizable variables.
    """
    def __init__(self, myclasswithoptimizablevariables):
        # doc needs to be initialized first
        # self.doc = doc

        self.ui_class = UIInterfaceClassWithOptimizableVariables(
                myclasswithoptimizablevariables)
        dict_from_query = self.ui_class.queryForDictionary()
        string_dict_from_query = self.ui_class.transformDictionaryForUI(
                dict_from_query)

        pprint(string_dict_from_query)

        filename = getRelativeFilePath(__file__, 'Qt/dlg_cwov_edit.ui')
        # this will create a Qt widget from our ui file
        self.form = FreeCADGui.PySideUic.loadUi(filename)
        self.readTableFromList(self.form.tableWidget_annotations,
                               [(k, v, True) for (k, v) in
                                string_dict_from_query["annotations"].items()])
        self.readTableFromList(self.form.tableWidget_variables,
                               string_dict_from_query["variables_list"])
        self.form.lineEdit_kind.setText(string_dict_from_query["kind"])
        self.form.lineEdit_name.setText(string_dict_from_query["name"])
        self.form.lineEdit_unique_id.setText(string_dict_from_query["unique_id"])

    def accept(self):
        dict_from_query = self.ui_class.queryForDictionary()
        string_dict_from_query = self.ui_class.transformDictionaryForUI(dict_from_query)

        var_list = self.writeTableToList(self.form.tableWidget_variables)
        annotations_dict = self.writeTableToList(
                self.form.tableWidget_annotations)
        string_dict_from_query["name"] = self.form.lineEdit_name.text()
        string_dict_from_query["variables_list"] = var_list
        string_dict_from_query["annotations"] = annotations_dict
        #self.ui_class.modifyFromDictionary(dict_from_query,
        #                                   transformation_dictionary)

        pprint(self.ui_class.queryForDictionary())

        FreeCADGui.Control.closeDialog()

    def reject(self):
        FreeCADGui.Control.closeDialog()

    def readTableFromList(self, mytable, mylist):
        """
        mylist contains triples of (name, value, modifyable)
        """
        mytable.clear()
        mytable.setRowCount(0)
        for (ind, (name, value, modifyable, var_type)) in enumerate(
                sorted(mylist, key=lambda x: x[0])):
            # sort list to get a reproducible table
            mytable.insertRow(ind)
            mytable.setItem(ind, 0, QTableWidgetItem(name))
            value_item = QTableWidgetItem(str(value))
            if not modifyable:
                value_item.setFlags(value_item.flags() & Qt.ItemIsEditable)
            mytable.setItem(ind, 1, value_item)

        header = mytable.horizontalHeader()

        try:
            # this is Qt4
            header.setResizeMode(0, QHeaderView.ResizeToContents)
            header.setResizeMode(1, QHeaderView.Stretch)
        except AttributeError:
            # this is Qt5
            header.setSectionResizeMode(0, QHeaderView.ResizeToContents)
            header.setSectionResizeMode(1, QHeaderView.Stretch)

    def writeTableToList(self, mytable):
        myvars = []
        for ind in range(mytable.rowCount()):
            var_name = mytable.item(ind, 0).text()
            var_value = float(mytable.item(ind, 1).text())
            var_modifyable = mytable.item(ind, 1).flags() != Qt.ItemIsEditable
            myvars.append((var_name, var_value, var_modifyable))
        return myvars
