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
#from PySide.QtGui import QTableWidgetItem

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
        print(dict_from_query)

        filename = getRelativeFilePath(__file__, 'Qt/dlg_cwov_edit.ui')
        # this will create a Qt widget from our ui file
        print(filename)
        self.form = FreeCADGui.PySideUic.loadUi(filename)

    def accept(self):
        FreeCADGui.Control.closeDialog()

    def reject(self):
        FreeCADGui.Control.closeDialog()

