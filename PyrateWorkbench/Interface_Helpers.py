# -*- coding: utf-8 -*-
"""
Created on Sat Nov 19 15:27:15 2016

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

from string import maketrans
import os

from Interface_Identifiers import *
from Interface_Checks import *

import FreeCAD, FreeCADGui

# general

def uuidToName(uuid):
    tab = maketrans('-','_')
    return str(uuid).lower().translate(tab)
    
def getRelativeFilePath(relativefilename, targetfile):
    return os.path.join(os.path.dirname(relativefilename), targetfile)
    
# collect optical system observers
    
def getAllOpticalSystemObservers(doc):
    return [obj for obj in doc.Objects if isOpticalSystemObserver(obj)]

# collect function objects
    
def getFunctionObjectsSubgroupFromOpticalSystemObserver(doc, os):
    fngroupname = os.NameFunctionsGroup
    fngroup = self.doc.getObject(fngroupname)
    return fngroup

def getFunctionObjectsFromOpticalSystemObserver(doc, os):
    subgroup = getFunctionObjectsSubgroupFromOpticalSystemObserver(doc, os)
    return subgroup.Group
    
# collect material catalogues

def getAllMaterialCatalogues(doc):
    return [obj for obj in doc.Objects if isMaterialCatalogue(obj)]
    
        
