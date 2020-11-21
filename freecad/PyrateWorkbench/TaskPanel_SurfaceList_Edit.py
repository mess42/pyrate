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

import FreeCADGui

from .Interface_Helpers import *
from .Interface_Checks import *

class SurfaceListTaskPanelEdit:
    def __init__(self, os):
        # doc needs to be initialized first
        self.os = os


        fn = getRelativeFilePath(__file__, 'Qt/dlg_surfacelist_edit.ui')


        # this will create a Qt widget from our ui file
        self.form = FreeCADGui.PySideUic.loadUi(fn)


    def accept(self):


        FreeCADGui.Control.closeDialog()
